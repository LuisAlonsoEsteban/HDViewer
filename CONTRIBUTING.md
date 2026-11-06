# Contributing to VisiumVisualizer

Thank you for helping to improve VisiumVisualizer. Bug reports, documentation fixes,
new features and new analysis methods are all welcome.

## Adding your method as a plugin

Segmentation and deconvolution methods are plugins: one folder with a manifest and a
script, with no changes to the viewer. Read [docs/plugins.md](docs/plugins.md); it
takes you from an empty folder to a method that runs in the interface.

To propose a method for the public instance at <https://hdviewer.unav.edu>:

1. Optionally, open a [*New method plugin* issue](https://github.com/LuisAlonsoEsteban/VisiumVisualizer/issues/new?template=method-plugin.yml)
   to discuss it first.
2. Fork the repository and add your plugin folder under
   `hd_viewer/startup/data/image/segmentation/api/plugins/` or
   `hd_viewer/analysis/deconvolution/api/plugins/`.
3. Go through the checklist in [docs/plugins.md](docs/plugins.md#before-you-share-checklist).
4. Open a pull request describing the method, its reference, licence and runtime
   (CPU/GPU, memory, time on a typical slide), and attach a screenshot of its result in
   the viewer.

After review and merge, the image workflow builds and publishes your plugin's image,
and the method becomes available to every user of the public instance, credited with
your citation.

## Reporting a problem

Open an [issue](https://github.com/LuisAlonsoEsteban/VisiumVisualizer/issues). Please
include what you did, what you expected, what happened, and the dataset (example name,
or the bin size and image size of your own data). Never attach private data.

## Development

```bash
pip install -r requirements.txt
uvicorn hd_viewer.api.main:app --port 8000      # http://localhost:8000
```

The architecture and the front-end workflow (Node.js 20+, `npm run dev`,
`npm run build`) are described under *For developers* in the [README](README.md).
Before opening a pull request, run the tests:

```bash
python -m unittest discover -s tests
cd frontend && npm test
```

By contributing, you agree that your contribution is released under the project's
[MIT License](LICENSE). A plugin may wrap a method under another licence, as long as
that licence allows redistribution and is stated with the plugin.
