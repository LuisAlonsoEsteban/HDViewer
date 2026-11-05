"""Barcodes in hover tooltips without sending one string per row.

Streamed rows carry no names. The tooltip formats the hovered row's position
(`$index`) with a CustomJSHover that looks the name up in a browser-side cache;
on a miss it shows "…" and asks the server (batched per animation tick) through
a small DataModel channel, then redraws the tooltip once the names arrive.

Row positions change only when rows are dropped or reloaded; the server then
bumps `version`, which empties the browser cache. Requests and answers carry the
version they refer to, so a name is never shown for the wrong row.
"""
import json

from bokeh.core.properties import Int, String
from bokeh.model import DataModel
from bokeh.models import CustomJS, CustomJSHover


class HoverNameChannel(DataModel):
    request = String(default="")
    response = String(default="")
    version = Int(default=0)


_FORMAT = """
const st = channel.__names ?? (channel.__names = {cache: new Map(), queued: new Set(), asked: new Set(),
                                                   version: null, seq: 0, timer: null})
if (st.version !== channel.version) {
    st.cache.clear(); st.queued.clear(); st.asked.clear(); st.version = channel.version
}
const name = st.cache.get(value)
if (name !== undefined) return name
if (!st.asked.has(value)) {
    st.queued.add(value)
    if (st.timer == null) {
        st.timer = requestAnimationFrame(() => {
            st.timer = null
            const rows = [...st.queued]
            st.queued.clear()
            rows.forEach((row) => st.asked.add(row))
            if (rows.length) channel.request = JSON.stringify({seq: ++st.seq, version: st.version, rows})
        })
    }
}
return "…"
"""

_RECEIVE = """
const st = channel.__names
if (st == null || channel.response === "") return
const message = JSON.parse(channel.response)
if (message.version !== st.version) return
for (const [row, name] of message.rows) {
    st.cache.set(row, name)
    st.asked.delete(row)
}
// Redraw the open tooltip at the last mouse position (no need to move the mouse)
for (const tool of tools) {
    const view = Bokeh.index.find_one(tool)
    if (view != null && view._current_sxy != null) view._inspect(...view._current_sxy)
}
"""


class HoverNames:
    """Server side of the lookup for one layer; `names()` gives the loaded rows'
    names in source order and `version()` the current row-position version."""

    FIELD = "$index"
    TEMPLATE = "$index{custom}"

    def __init__(self, names, version):
        self._names, self._version = names, version
        self.channel = HoverNameChannel()
        self.formatter = CustomJSHover(args=dict(channel=self.channel), code=_FORMAT)
        self._receive = CustomJS(args=dict(channel=self.channel, tools=[]), code=_RECEIVE)
        self.channel.js_on_change("response", self._receive)
        self.channel.on_change("request", self._answer)

    @property
    def formatters(self):
        return {self.FIELD: self.formatter}

    def bind(self, tool):
        """Redraw `tool`'s tooltip when names arrive."""
        if tool not in self._receive.args["tools"]:
            self._receive.args = dict(self._receive.args, tools=[*self._receive.args["tools"], tool])

    def sync_version(self):
        version = self._version()
        if self.channel.version != version:
            self.channel.version = version

    def lookup(self, version, rows):
        """[(row, name)] for `rows` of the given version (none if it is stale)."""
        if version != self._version():
            return []
        names = self._names()
        return [(int(row), str(names[row])) for row in rows
                if isinstance(row, int) and 0 <= row < len(names)]

    def _answer(self, attr, old, new):
        if not new:
            return
        try:
            message = json.loads(new)
            rows = self.lookup(message.get("version"), message.get("rows", [])[:1000])
        except (ValueError, TypeError, AttributeError):
            return
        self.channel.response = json.dumps(dict(seq=message.get("seq"), version=self._version(), rows=rows))
