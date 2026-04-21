
from . import Panel

from bokeh.models import Div

class MetadataPanel(Panel):

    def __init__(self):
        super().__init__(name="metadata")

    @staticmethod
    def parse_metadata(geom):
        return f"""
        <div>
            <p><strong>{geom.alias}s:</strong><br>
            Total {geom.alias.lower()}s: {geom.sample.n_geometries}<br>
            Not Sequenced: {geom.sample.n_not_sequenced}</p>
        </div>
        """

    def populate(self, geometries):
        text = ""
        for geom in geometries:
            text += MetadataPanel.parse_metadata(geom)
        
        self.set_children([Div(text=f"<div><h3>Metadata</h3>{text}</div>")])



        