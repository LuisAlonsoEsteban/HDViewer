/* Landing, setup, and workspace share one live Bokeh session. */
(() => {
  "use strict";
  const byId = (id) => document.getElementById(id);
  const screens = {landing: "landing-container", setup: "startup-container", workspace: "main-container"};
  let workspaceLoaded = false;
  let splits = [];

  function resizePlots() {
    requestAnimationFrame(() => window.dispatchEvent(new Event("resize")));
  }

  function initializeSplits() {
    if (splits.length || typeof window.Split !== "function" || window.innerWidth <= 700) return;
    const settings = JSON.parse(byId("layout-settings").textContent);
    const common = {gutterSize: 7, onDragEnd: resizePlots};
    splits = [
      Split(["#left", "#center", "#right"], {...common, ...settings.columns}),
      Split(["#metadata-container", "#widgets-container"], {...common, ...settings.metadata, direction: "vertical"}),
      Split(["#map-container", "#de-results-container"], {...common, ...settings.map, direction: "vertical"}),
      Split(["#top-right", "#bottom-right"], {...common, ...settings.right, direction: "vertical"}),
    ];
  }

  function showScreen(screen) {
    for (const [name, id] of Object.entries(screens)) byId(id).hidden = name !== screen;
    document.body.dataset.screen = screen;
    document.querySelector('.skip-link').href = screen === 'workspace' ? '#center' : screen === 'setup' ? '#startup-container' : '#landing-main';
    if (screen === "workspace") initializeSplits();
    window.scrollTo(0, 0);
    resizePlots();
  }

  window.showMainApp = () => {
    workspaceLoaded = true;
    document.querySelectorAll("[data-launch-label]").forEach((label) => label.textContent = "Resume workspace");
    showScreen("workspace");
  };

  document.addEventListener("click", (event) => {
    const action = event.target.closest("[data-action]")?.dataset.action;
    if (action === "launch") showScreen(workspaceLoaded ? "workspace" : "setup");
    if (action === "home") showScreen("landing");

    const accordion = event.target.closest(".widget-title");
    if (accordion) {
      const expanded = accordion.getAttribute("aria-expanded") !== "true";
      accordion.setAttribute("aria-expanded", String(expanded));
      byId(accordion.getAttribute("aria-controls")).classList.toggle("visible", expanded);
      resizePlots();
    }

    const preview = event.target.closest("[data-preview]");
    if (preview) {
      document.querySelectorAll("[data-preview]").forEach((button) => button.setAttribute("aria-pressed", String(button === preview)));
      byId("tissue-preview").dataset.view = preview.dataset.preview;
      const captions = {tissue: "Original tissue image", bins: "Illustrative spatial grid", regions: "Illustrative regions of interest"};
      byId("preview-caption").textContent = captions[preview.dataset.preview];
    }
  });

  byId("copy-citation").addEventListener("click", async () => {
    const citation = byId("software-citation");
    const status = byId("citation-status");
    try {
      await navigator.clipboard.writeText(citation.textContent.trim());
      status.textContent = "Citation copied";
    } catch (_) {
      const selection = window.getSelection();
      const range = document.createRange();
      range.selectNodeContents(citation);
      selection.removeAllRanges();
      selection.addRange(range);
      status.textContent = "Citation selected — copy with Ctrl/Cmd+C";
    }
  });

  const mobile = window.matchMedia("(max-width: 700px)");
  mobile.addEventListener("change", () => {
    splits.forEach((split) => split.destroy());
    splits = [];
    if (!byId(screens.workspace).hidden) initializeSplits();
    resizePlots();
  });
})();
