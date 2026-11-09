// The public landing page (ported from the Bokeh interface's templates/index.html).
import { useState } from "react";
import "../landing/landing.css";
import preview from "../landing/tissue-preview.svg?raw";
import { useApp } from "../store/app";

// Files of hd_viewer/static, served at static/assets/ (relative: works under any prefix)
const STATIC = "static/assets";

export function BrandMark({ className = "brand-mark" }: { className?: string }) {
  return (
    <svg className={className} viewBox="0 0 36 36" fill="none" aria-hidden="true">
      <path d="M18 2 32 10v16L18 34 4 26V10Z" stroke="currentColor" strokeWidth="1.5" />
      <path d="m11 14 7-4 7 4v8l-7 4-7-4Z" fill="currentColor" opacity=".14" />
      <circle cx="18" cy="10" r="2.4" fill="currentColor" />
      <circle cx="11" cy="22" r="2.4" fill="currentColor" />
      <circle cx="25" cy="22" r="2.4" fill="currentColor" />
      <path d="m18 12-6 8m7-8 5 8m-10 2h8" stroke="currentColor" strokeWidth="1.3" />
    </svg>
  );
}

const PREVIEWS = {
  tissue: "Original tissue image",
  bins: "Illustrative spatial grid",
  regions: "Illustrative regions of interest",
} as const;
type Preview = keyof typeof PREVIEWS;

const external = { target: "_blank", rel: "noopener noreferrer" };

export function Landing() {
  const { config, workspace, setScreen } = useApp();
  const branding = config?.branding ?? {};
  const tutorial = config?.tutorial_url ?? `${STATIC}/tutorial.html`;
  const useCase = config?.use_case_url ?? `${STATIC}/biological-use-case.html`;
  const [view, setView] = useState<Preview>("tissue");
  const [citationStatus, setCitationStatus] = useState("");

  // "Launch" resumes the loaded workspace, as in the old interface
  const launch = () => setScreen(workspace ? "workspace" : "loading");
  const launchLabel = (label: string) => (workspace ? "Resume workspace" : label);

  const copyCitation = async () => {
    try {
      await navigator.clipboard.writeText(branding.citation ?? "");
      setCitationStatus("Citation copied");
    } catch {
      const node = document.getElementById("software-citation");
      if (node) {
        const range = document.createRange();
        range.selectNodeContents(node);
        window.getSelection()?.removeAllRanges();
        window.getSelection()?.addRange(range);
      }
      setCitationStatus("Citation selected — copy with Ctrl/Cmd+C");
    }
  };

  return (
    <div className="landing">
      <a className="skip-link" href="#landing-main">Skip to content</a>
      <header className="site-header wrap">
        <a className="brand" href="#" aria-label={`${branding.name ?? "VisiumVisualizer"} home`}>
          <BrandMark />{branding.name ?? "VisiumVisualizer"}
        </a>
        <nav className="site-nav" aria-label="Main navigation">
          <a className="desktop-link" href="#features">Explore</a>
          <a className="desktop-link" href="#workflow">Workflow</a>
          <a className="desktop-link" href="#plugins">Plugins</a>
          <a className="desktop-link" href="#research">Cite &amp; contact</a>
          {branding.repository && <a href={branding.repository} {...external}>GitHub ↗</a>}
          <button className="button small" onClick={launch}>{launchLabel("Open workspace")} <span aria-hidden="true">↗</span></button>
        </nav>
      </header>

      <main id="landing-main" tabIndex={-1}>
        <section className="hero wrap" aria-labelledby="hero-title">
          <div>
            <div className="eyebrow">Spatial transcriptomics · Visium HD</div>
            <h1 id="hero-title">See the tissue.<br />Explore the <em>biology.</em></h1>
            <p className="intro">
              Bring gene expression, cell populations, and tissue neighborhoods into one view. An interactive
              workspace for asking spatial questions, with the original tissue image as your guide.
            </p>
            <div className="hero-actions">
              <button className="button" onClick={launch}>{launchLabel("Explore your data")} <span aria-hidden="true">↗</span></button>
              <a className="button secondary" href={tutorial} {...external}>Read the guide <span aria-hidden="true">↗</span></a>
            </div>
            <p className="hero-note">
              Bring your Visium HD files, or start with a configured example. Free and open to everyone, no login
              required; uploaded data stays private to your session and is deleted when it ends.
            </p>
          </div>
          <figure className="tissue-figure">
            <div className="preview-window">
              <div className="preview-toolbar">
                <span><span className="dots" aria-hidden="true"><i /><i /><i /></span> HUMAN COLORECTAL TISSUE</span>
                <span>P1CRC / H&amp;E</span>
              </div>
              <div className="preview-stage" data-view={view}>
                <img src={`${STATIC}/assets/p1crc-tissue.png`} width={600} height={496}
                     alt="Hematoxylin and eosin tissue image from the human colorectal cancer example" />
                {/* Static illustrative overlay (bins and regions), our own markup */}
                <div dangerouslySetInnerHTML={{ __html: preview }} />
                <div className="preview-tools" aria-hidden="true"><span>+</span><span>−</span><span>⌖</span></div>
                <span className="preview-scale">{PREVIEWS[view]}</span>
              </div>
              <div className="preview-tabs" role="group" aria-label="Illustrative tissue preview">
                {(Object.keys(PREVIEWS) as Preview[]).map((key) => (
                  <button key={key} aria-pressed={view === key} onClick={() => setView(key)}>
                    {{ tissue: "Tissue", bins: "Spatial bins", regions: "Regions" }[key]}
                  </button>
                ))}
              </div>
            </div>
            <figcaption>
              <span>One tissue. Multiple ways to explore.</span>
              <a href={useCase} {...external}>View the biological example ↗</a>
            </figcaption>
          </figure>
        </section>

        <section className="affiliations wrap" aria-label="Research affiliations">
          <p>Connecting computational biology<br />and biomedical research.</p>
          <div className="affiliation-logos">
            <a href="https://www.unav.edu/" {...external}>
              <img src={`${STATIC}/assets/university-navarra.svg`} alt="University of Navarra" width={205} height={43} loading="lazy" />
            </a>
            <a href="https://tecnun.unav.edu/" {...external}>
              <img src={`${STATIC}/assets/tecnun.jpg`} alt="Tecnun Universidad de Navarra — School of Engineering" width={540} height={360} loading="lazy" />
            </a>
            <a href="https://cima.cun.es/" {...external}>
              <img src={`${STATIC}/assets/cima.jpg`} alt="Cima Universidad de Navarra" width={738} height={273} loading="lazy" />
            </a>
          </div>
        </section>

        <section id="features" className="section wrap" aria-labelledby="features-title">
          <div className="section-heading">
            <div>
              <div className="eyebrow">Explore in context</div>
              <h2 id="features-title">From spatial bins to biological questions.</h2>
            </div>
            <p>Keep the image, the expression, and the populations connected as you move between scales.</p>
          </div>
          <div className="feature-list">
            <article className="feature">
              <span className="number">01 / TISSUE &amp; EXPRESSION</span>
              <h3>Find where a gene belongs.</h3>
              <p>Browse gene expression on the tissue. Move between native bins, aggregated views, pseudospots, and segmented cell layers.</p>
              <a href={`${tutorial}#app-color`} {...external}>Explore the map ↗</a>
            </article>
            <article className="feature">
              <span className="number">02 / CELLS &amp; POPULATIONS</span>
              <h3>Give each region context.</h3>
              <p>Compare segmentation methods, cluster expression profiles, and estimate cell-type composition with deconvolution.</p>
              <a href={`${tutorial}#app-clustering`} {...external}>Understand the populations ↗</a>
            </article>
            <article className="feature">
              <span className="number">03 / REGIONS &amp; COMPARISONS</span>
              <h3>Ask a spatial question.</h3>
              <p>Draw a region, start from a barcode, or build a neighborhood around a group. Compare expression and composition across regions.</p>
              <a href={useCase} {...external}>Follow a biological use case ↗</a>
            </article>
          </div>
        </section>

        <section id="workflow" className="workflow-section section" aria-labelledby="workflow-title">
          <div className="wrap workflow">
            <div>
              <div className="eyebrow">A connected workflow</div>
              <h2 id="workflow-title">Start with your tissue.<br />Follow your question.</h2>
              <p>Use your own spatial data or explore an example. Each layer stays connected to the tissue, so your comparisons keep their spatial context.</p>
              <button className="button small" onClick={launch}>{launchLabel("Set up a dataset")} ↗</button>
            </div>
            <ol className="workflow-steps">
              <li><span className="step-marker">01</span><div><h3>Load &amp; prepare</h3><p>Choose an expression matrix, tissue positions, and an image. Add a single-cell reference and select optional preprocessing.</p></div></li>
              <li><span className="step-marker">02</span><div><h3>Explore &amp; define</h3><p>Find genes and cell populations. Use linked views, native barcodes, and regions of interest to define your comparison.</p></div></li>
              <li><span className="step-marker">03</span><div><h3>Compare &amp; interpret</h3><p>Examine markers, differential expression, and population composition. Follow the worked colorectal cancer example for a complete analysis.</p></div></li>
            </ol>
          </div>
        </section>

        <section id="plugins" className="section wrap" aria-labelledby="plugins-title">
          <div className="section-heading">
            <div>
              <div className="eyebrow">Open to new methods</div>
              <h2 id="plugins-title">Bring your own method.</h2>
            </div>
            <p>
              Every segmentation and deconvolution method here is a plugin. Add yours the same way, with no
              changes to the viewer.
            </p>
          </div>
          <div className="feature-list">
            <article className="feature">
              <span className="number">ONE FOLDER</span>
              <h3>A manifest and a script.</h3>
              <p>
                <code>manifest.yaml</code> describes your method and its settings; the interface builds the form from
                it. <code>run.py</code> reads the inputs and writes the result, with a small helper library.
              </p>
            </article>
            <article className="feature">
              <span className="number">ANY STACK</span>
              <h3>Your dependencies, isolated.</h3>
              <p>
                Plugins run in their own container (Docker or Apptainer) or Python environment: Python, R or GPU
                code alongside the bundled methods, without conflicts.
              </p>
            </article>
            <article className="feature">
              <span className="number">SAME RESULTS VIEW</span>
              <h3>Compared side by side.</h3>
              <p>
                A segmentation plugin adds its own cells layer; a deconvolution plugin adds runs to the
                Deconv. tab. Both can be compared with the bundled methods, and your citation is shown with them.
              </p>
            </article>
          </div>
          <div className="hero-actions">
            <a className="button small" href={branding.plugins_guide} {...external}>Read the plugin guide ↗</a>
            <a className="button secondary small" href={branding.plugin_proposal} {...external}>Propose a method ↗</a>
            <a className="button secondary small" href={`${tutorial}#app-plugins`} {...external}>How plugins appear in the app ↗</a>
          </div>
        </section>

        <section id="research" className="section wrap" aria-labelledby="research-title">
          <div className="research-grid">
            <div>
              <div className="eyebrow">Research &amp; reproducibility</div>
              <h2 id="research-title">Built to be used.<br />Easy to reference.</h2>
              <div className="citation-box">
                <p id="software-citation">{branding.citation}</p>
                <div className="citation-actions">
                  <button className="button secondary small" onClick={copyCitation}>Copy citation</button>
                  <a href={`${STATIC}/citation.cff`} download>Download citation ↧</a>
                  <span role="status" aria-live="polite">{citationStatus}</span>
                </div>
              </div>
              <p className="source-note">
                Please also cite the methods and datasets used in your analysis. The colorectal example follows{" "}
                <a href={branding.data_paper} {...external}>Oliveira et al., Nature Genetics (2025)</a>; this is the
                source-data paper, separate from the software reference above.
              </p>
            </div>
            <aside className="contact-panel" id="contact" aria-labelledby="contact-title">
              <div className="eyebrow">Let’s connect</div>
              <h3 id="contact-title">Questions, ideas, or a useful finding?</h3>
              <p>Share a reproducible issue, suggest an improvement, or get in touch about the project.</p>
              <p>
                <span className="contact-name">{branding.author}</span><br />
                <a href={branding.affiliation_url} {...external}>{branding.affiliation}</a><br />
                <a href={`mailto:${branding.email}`}>{branding.email}</a>
              </p>
              <div className="hero-actions">
                <a className="button small" href={branding.issues} {...external}>GitHub issues ↗</a>
                <a className="button secondary small" href={branding.repository} {...external}>View source ↗</a>
              </div>
            </aside>
          </div>
        </section>
      </main>

      <footer className="closing wrap">
        <span>{branding.name ?? "VisiumVisualizer"} · Spatial biology, in context.</span>
        <nav aria-label="Footer">
          <a href={tutorial} {...external}>Documentation</a>
          <a href={branding.plugins_guide} {...external}>Plugins</a>
          <a href={branding.repository} {...external}>GitHub ↗</a>
          <a href={`mailto:${branding.email}`}>Contact</a>
        </nav>
      </footer>
    </div>
  );
}
