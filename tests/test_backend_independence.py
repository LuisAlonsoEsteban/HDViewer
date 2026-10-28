"""Backend dependency boundary and application-to-backend configuration contracts."""
import ast
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import unittest
from unittest.mock import Mock, patch

import numpy as np

from bokeh_server import backend_config, settings
from bokeh_server.components.image import ImageComponent
from bokeh_server.components.utils.selection import native_selection, selection_indices
from src.data.image.image import ImageProcessor
from src.data.image import transform as spatial
from src.operations.bins2spots import Bins2Spots
from src.sample import Bins, Cells
from src.segmentation.service import SegmentationService
from src.segmentation.runners import make_runner
from src.utils import comparison
from tests.test_spatial_alignment import fixture

ROOT = Path(__file__).resolve().parents[1]


class BackendBoundaryTests(unittest.TestCase):
    def test_src_has_no_frontend_imports(self):
        for path in (ROOT/'src').rglob('*.py'):
            for node in ast.walk(ast.parse(path.read_text())):
                names = [a.name for a in node.names] if isinstance(node,ast.Import) else [node.module or ''] if isinstance(node,ast.ImportFrom) else []
                self.assertFalse(any(n=='bokeh_server' or n.startswith('bokeh_server.') for n in names), str(path))

    def test_all_backend_modules_and_analysis_work_when_frontend_imports_are_blocked(self):
        code = textwrap.dedent('''
            import importlib, importlib.abc, pkgutil, sys
            class NoFrontend(importlib.abc.MetaPathFinder):
                def find_spec(self, fullname, path=None, target=None):
                    if fullname == 'bokeh_server' or fullname.startswith('bokeh_server.'):
                        raise AssertionError('Backend attempted frontend import: ' + fullname)
            sys.meta_path.insert(0, NoFrontend())
            import src
            for info in pkgutil.walk_packages(src.__path__, 'src.'):
                importlib.import_module(info.name)
            import anndata as ad
            import numpy as np
            import pandas as pd
            from src.utils.deconvolution import deconvolve
            from src.sample._sample import Sample
            reference = ad.AnnData(np.array([[10.,0.],[0.,10.]]),
                obs=pd.DataFrame({'cell_type':['a','b']},index=['r1','r2']),var=pd.DataFrame(index=['g1','g2']))
            target = ad.AnnData(np.array([[6.,4.],[3.,7.]]),
                obs=pd.DataFrame(index=['s1','s2']),var=reference.var.copy())
            result, genes = deconvolve(target, reference, 'cell_type')
            np.testing.assert_allclose(result[['a','b']].sum(axis=1), 1.)
            sample = Sample(pd.DataFrame({'index':['s1','s2']}), target)
            assert sample.cluster_col == 'leiden'
            assert not any(n=='bokeh_server' or n.startswith('bokeh_server.') for n in sys.modules)
        ''')
        result = subprocess.run([sys.executable,'-c',code],cwd=ROOT,capture_output=True,text=True,timeout=120)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)


class BackendConfigurationTests(unittest.TestCase):
    def test_nondefault_sample_settings_and_aggregate_inheritance(self):
        _, transform, _, original = fixture()
        with patch.multiple(settings,GENE_CACHE_SIZE=2,MIN_GENE_CELLS=1,NORMALIZATION_TARGET_SUM=123.,
                            CLUSTER_COL='types',CLUSTER_PREFIX='types_',AGGREGATED_PREFIX='level:',
                            AGGREGATE_OBS_PREFIX='block:',AGGREGATION_FACTORS=(1,2),MAX_CACHE_ENTRIES=2):
            sample = Bins(original.gdf,original.adata,spatial_transform=transform,pixels_per_bin=4,
                          **backend_config.sample_kwargs('bins'))
        self.assertEqual(sample.get_resolution_level_for_factor(2),'level:2')
        self.assertEqual(sample.get_cluster_colname(2),'types_2')
        names = sample.aggregation_names(2)
        first = names.iloc[0]
        self.assertTrue(first.startswith('block:'))
        selected = native_selection(sample,{first})
        self.assertEqual(selected,set(names.index[names==first]))
        self.assertEqual(selection_indices(sample,selected,[first]),[0])
        aggregate = sample._get_aggregated_sample(2)
        self.assertEqual(aggregate._sample_options,sample._sample_options)
        self.assertEqual(aggregate.cluster_col,'types')
        self.assertEqual(aggregate.GENE_CACHE_SIZE,2)
        self.assertEqual(aggregate.normalization_target_sum,123.)
        with patch('scanpy.pp.filter_genes') as filt, patch('scanpy.pp.normalize_total') as norm, \
             patch('scanpy.pp.log1p'), patch('scanpy.pp.highly_variable_genes') as hvg:
            hvg.side_effect=lambda a: a.var.__setitem__('highly_variable',True)
            aggregate.prepare_adata_processed()
            self.assertEqual(filt.call_args.kwargs,dict(min_cells=1))
            self.assertEqual(norm.call_args.kwargs,dict(target_sum=123.))
            self.assertIsNot(norm.call_args.args[0],aggregate.adata)
        with patch.object(aggregate,'compute_embedding',return_value=np.zeros((aggregate.adata.n_obs,2))) as embed:
            data = sample.get_embedding_dataframe('pca',resolution_level='level:2',n_pcs=7,n_neighbors=4)
            embed.assert_called_once_with('pca',n_pcs=7,n_neighbors=4)
            self.assertEqual(len(data),aggregate.adata.n_obs)
        with patch.multiple(settings,LOD_PREFIX='detail:',LOD_BASE_TOLERANCE=3.,LOD_MAX_LEVEL=6,LOD_TOLERANCE_PX=2.,LOD_HYSTERESIS=.2):
            cells=Cells(original.gdf,original.adata,**backend_config.sample_kwargs('cells'))
        self.assertEqual(cells.get_tolerance('detail:2'),12.)
        self.assertEqual(Cells.LOD_PREFIX,'lod:')

    def test_application_startup_passes_configured_backend_arguments(self):
        from tests.test_spatial_alignment import AlignedIntegrationTests
        from src.sample import Bins as BackendBins
        with patch.multiple(settings,TILE_SIZE=32,GRID_FIT_MAX_POINTS=17,
                            SPOT_BUFFER_QUAD_SEGS=8,GENE_CACHE_SIZE=3,MIN_GENE_CELLS=2), \
             patch('src.data.image.image.ImageProcessor',wraps=ImageProcessor) as image_factory, \
             patch('src.operations.bins2spots.Bins2Spots',wraps=Bins2Spots) as spots_factory, \
             patch('bokeh_server.panels.app.app.Bins',wraps=BackendBins) as bins_factory:
            AlignedIntegrationTests('test_startup_builds_aligned_image_before_native_layers').test_startup_builds_aligned_image_before_native_layers()
            self.assertEqual(image_factory.call_args.kwargs['tile_size'],32)
            self.assertEqual(spots_factory.call_args.kwargs['spot_buffer_quad_segs'],8)
            self.assertEqual(spots_factory.call_args.kwargs['grid_fit_max_points'],17)
            self.assertEqual(bins_factory.call_args.kwargs['gene_cache_size'],3)
            self.assertEqual(bins_factory.call_args.kwargs['min_gene_cells'],2)

    def test_grid_converter_tile_and_server_options(self):
        positions, _, _, sample = fixture()
        with patch.multiple(settings,GRID_FIT_MAX_POINTS=17,GRID_MAX_CONDITION_NUMBER=1e9,
                            SPOT_BUFFER_QUAD_SEGS=8,TILE_SIZE=64):
            with patch.object(spatial,'fit_grid_basis',wraps=spatial.fit_grid_basis) as fit:
                transform=spatial.SpatialTransform.fit(positions,(140,180,3),**backend_config.grid_fit_kwargs())
                self.assertEqual(fit.call_args.kwargs,dict(max_points=17,max_condition_number=1e9))
            b2s=Bins2Spots(sample.adata,positions,bins_size=2,**backend_config.bins_to_spots_kwargs())
            self.assertEqual(b2s.converter.visiumHD.grid_fit_max_points,17)
            self.assertEqual(b2s.spot_buffer_quad_segs,8)
            image=ImageProcessor(None,'unused',spatial_transform=transform,**backend_config.image_kwargs())
            self.assertEqual(image.tile_size,64)
            self.assertIsNone(image._image)
        processor=Mock()
        processor.deploy_server.return_value='key'
        with patch.object(settings,'SERVER_HOST','127.0.0.2'):
            ImageComponent(processor,port=9876).deploy_server()
        processor.deploy_server.assert_called_once_with(host='127.0.0.2',port=9876)

    def test_segmentation_paths_deep_merge_and_environment_override(self):
        with tempfile.TemporaryDirectory() as root:
            folder=Path(root)
            config=folder/'config.yaml'
            config.write_text('runner: local\nthreads_per_run: 2\ndocker:\n  binary: special-docker\napptainer:\n  images:\n    watershed: example.sif\n')
            plugins=folder/'plugins'; plugins.mkdir()
            with patch.multiple(settings,SEGMENTATION_CONFIG_PATH=config,PLUGINS_DIR=plugins,
                                SDK_DIR=folder/'sdk',SEGMENTATION_ROOT=folder,DEFAULT_THREADS=7,
                                SEGMENTATION_MAX_PARALLEL=3,SEGMENTATION_GPU=True), \
                 patch.dict(os.environ,{'VHD_SEG_RUNNER':''}):
                options=backend_config.segmentation_kwargs()
                service=SegmentationService(**options)
                self.assertEqual(service.plugins,{})
                self.assertEqual(service.runner.name,'local')
                self.assertEqual(service.runner.threads,2)
                self.assertTrue(service.runner.gpu)
                self.assertEqual(service._slots._value,3)
                self.assertEqual(service.runner.sdk_dir,folder/'sdk')
                self.assertEqual(options['config']['docker']['binary'],'special-docker')
                self.assertEqual(options['config']['apptainer']['images_dir'],settings.SEGMENTATION_IMAGES_DIR)
                self.assertEqual(options['config']['apptainer']['images'],{'watershed':'example.sif'})
                with patch.dict(os.environ,{'VHD_SEG_RUNNER':'docker'}):
                    runner=make_runner(options['config'],sdk_dir=options['sdk_dir'],segmentation_root=folder)
                    self.assertEqual(runner.name,'docker')
                    self.assertEqual(runner.segmentation_root,folder)

    def test_application_analysis_calls_pass_configured_values(self):
        from bokeh.events import ButtonClick
        from bokeh.models import Button
        from bokeh_server.components.utils.loading import LoadingLogger
        from bokeh_server.panels.app.embedding import EmbeddingPanel
        from tests.test_deconvolution import reference, ImmediateDoc
        from tests.test_param_form import make_layer
        component, _ = make_layer()
        def work_only(logger, title, work, *args, **kwargs):
            return work()
        with patch.multiple(settings,N_PCS=7,N_NEIGHBORS=4,DECONV_BATCH_SIZE=2), \
             patch.object(LoadingLogger,'run',work_only):
            with patch.object(component.sample,'run_clustering',return_value='leiden') as clustering:
                component.clustering_manager.compute(.65)
                self.assertEqual(clustering.call_args.kwargs['n_pcs'],7)
                self.assertEqual(clustering.call_args.kwargs['n_neighbors'],4)
                self.assertEqual(clustering.call_args.kwargs['resolution'],.65)
            panel=EmbeddingPanel(ImmediateDoc())
            panel.set_components([component])
            with patch.object(component.sample,'get_embedding_dataframe') as embedding:
                panel._run_embedding(component.key)
                self.assertEqual(embedding.call_args.kwargs['n_pcs'],7)
                self.assertEqual(embedding.call_args.kwargs['n_neighbors'],4)
            component.deconvolution_manager.reference=reference()
            widgets=component.deconvolution_manager.get_widgets(ImmediateDoc())
            run=next(w for w in widgets if isinstance(w,Button) and w.label=='Run deconvolution')
            with patch('bokeh_server.components.utils.deconvolution.deconvolve') as deconvolution:
                run._trigger_event(ButtonClick(run))
                self.assertEqual(deconvolution.call_args.kwargs['batch_size'],2)

    def test_comparison_options_reach_private_calculations(self):
        _, _, _, sample = fixture()
        with patch.object(comparison,'_log_normalize',wraps=comparison._log_normalize) as norm, \
             patch.object(comparison,'_profile_metrics',wraps=comparison._profile_metrics) as profile, \
             patch.object(comparison,'_knn_weights',wraps=comparison._knn_weights) as knn, \
             patch.object(comparison,'_tissue_hull',wraps=comparison._tissue_hull) as hull, \
             patch.object(comparison,'_morans_i',wraps=comparison._morans_i) as morans:
            result=comparison.compare_samples(sample,sample,genes=['gene','one'],n_top_svg=1,
                normalization_target_sum=123.,chunk_size=3,k_neighbors=2,hull_ratio=.4)
            self.assertIsNotNone(result['local'])
            self.assertEqual(norm.call_args.args[1],123.)
            self.assertEqual(profile.call_args.kwargs['chunk_size'],3)
            self.assertEqual(knn.call_args.kwargs['k'],2)
            self.assertEqual(hull.call_args.kwargs['ratio'],.4)
            self.assertTrue(all(c.kwargs.get('chunk_size')==3 for c in morans.call_args_list))


if __name__ == '__main__':
    unittest.main()
