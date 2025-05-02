
from functools import partial

class Widget:

    def __init__(self, component, callback, item):
        self.component = component
        self.callback = callback
        self.item = item

    def connect(self, source=None):
        self.component.on_change(self.item, partial(self.callback, source=source))


