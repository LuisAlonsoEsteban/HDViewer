
"""
var updateData = function () {
    var xhr = new XMLHttpRequest();
    var left = xr.start;
    var right = xr.end;
    var top = yr.start;
    var bottom = yr.end;
    var url = `/get_data?left=${left}&right=${right}&top=${top}&bottom=${bottom}`;
    xhr.open('GET', url);
    xhr.onload = function() {
        if (xhr.status === 200) {
            var newData = JSON.parse(xhr.responseText);
            source.data = newData;
            rg.visible = true;  // Ensure that the rectangle glyphs are visible
        }
    };
    xhr.send();
}
"""

callback = CustomJS(args=dict(rg=rect_glyph, xr=p.x_range, yr=p.y_range, source=source), code=updateData.toString() + " updateData();");

p.x_range.js_on_change('end', callback);
p.y_range.js_on_change('end', callback);
