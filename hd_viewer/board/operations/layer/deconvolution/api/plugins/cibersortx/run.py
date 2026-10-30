"""CIBERSORTx uses its official container, orchestrated by the host service.

The host converts AnnData to the official text format and executes the image
through the same Docker/Apptainer runner used by other plugins. This avoids
requiring a container runtime inside a container.
"""
if __name__ == '__main__':
    raise SystemExit('Use DeconvolutionService with Docker or Apptainer for CIBERSORTx.')
