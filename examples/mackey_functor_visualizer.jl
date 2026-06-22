### A Pluto.jl notebook ###
# v0.20.25

using Markdown
using InteractiveUtils

# ╔═╡ 11111111-1111-1111-1111-111111111111
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))

    using GAP
    using AbstractAlgebra
    using MackeyFunctors
end

# ╔═╡ 22222222-2222-2222-2222-222222222222
md"""
# Mackey Functor Visualizer

Edit the `M = ...` cell to visualize a different validated `MackeyFunctor`.
The diagram shows subgroup names and values. Click nodes for subgroup metadata;
click two comparable nodes to inspect composite restriction and transfer maps.
"""

# ╔═╡ 33333333-3333-3333-3333-333333333333
begin
    C2 = GAP.Globals.CyclicGroup(2)
    subs = Vector{GapObj}(GAP.Globals.AllSubgroups(C2))
    e = subs[1]
    g = GAP.Globals.GeneratorsOfGroup(C2)[1]

    values = IdDict{GapObj, Vector{Int64}}(
        e => [0],
        C2 => [0, 0],
    )

    restrictions = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (e, C2) => [1 2],
    )

    transfers = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (e, C2) => [0; 1;;],
    )

    conjugations = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (e, g) => [1;;],
    )

    M = MackeyFunctor(C2, values, restrictions, transfers, conjugations)
end

# ╔═╡ 44444444-4444-4444-4444-444444444444
function mackey_visualizer_html(M::MackeyFunctor)
    data_json = visualizer_json(M)

    return HTML("""
<div id="mackey-visualizer" style="font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; color: #000; font-size: 16px;">
  <style>
    #mackey-visualizer .mv-shell {
      position: relative;
      width: 100%;
      min-width: 0;
    }
    #mackey-visualizer .mv-toolbar {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      align-items: center;
      justify-content: space-between;
      margin-bottom: 12px;
    }
    #mackey-visualizer .mv-toolbar-group {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }
    #mackey-visualizer button {
      border: 1px solid #000;
      background: #fff;
      color: #000;
      min-height: 32px;
      border-radius: 999px;
      padding: 0 14px;
      font: inherit;
      font-size: 15px;
      cursor: pointer;
    }
    #mackey-visualizer button:hover {
      box-shadow: inset 0 0 0 1px #000;
    }
    #mackey-visualizer button.active {
      background: #000;
      color: #fff;
    }
    #mackey-visualizer .mv-graph {
      width: 100%;
      min-width: 0;
      min-height: min(760px, 72vh);
      overflow: auto;
      border: 1px solid #000;
      background: #fff;
      box-sizing: border-box;
    }
    #mackey-visualizer.expanded {
      position: fixed;
      inset: 16px;
      z-index: 1000;
      padding: 16px;
      background: #fff;
      border: 1px solid #000;
      box-shadow: 0 24px 64px rgba(0, 0, 0, 0.22);
      box-sizing: border-box;
    }
    #mackey-visualizer.expanded .mv-shell {
      height: calc(100vh - 96px);
    }
    #mackey-visualizer.expanded .mv-graph {
      height: 100%;
      min-height: 0;
    }
    #mackey-visualizer svg {
      display: block;
      max-width: none;
      background: #fff;
    }
    #mackey-visualizer .mv-panel {
      position: absolute;
      top: 14px;
      right: 14px;
      width: min(440px, calc(100% - 28px));
      max-height: calc(100% - 28px);
      overflow: auto;
      border: 1px solid #000;
      padding: 16px 18px;
      background: #fff;
      box-shadow: 0 16px 36px rgba(0, 0, 0, 0.18);
      opacity: 0;
      transform: translateY(-8px);
      pointer-events: none;
      transition: opacity 120ms ease, transform 120ms ease;
      box-sizing: border-box;
      z-index: 5;
    }
    #mackey-visualizer .mv-panel.open {
      opacity: 1;
      transform: translateY(0);
      pointer-events: auto;
    }
    #mackey-visualizer .mv-panel h3 {
      margin: 0 0 14px;
      font-size: 20px;
      line-height: 1.15;
    }
    #mackey-visualizer .mv-panel h4 {
      margin: 18px 0 8px;
      font-size: 16px;
      line-height: 1.2;
    }
    #mackey-visualizer .mv-panel dl {
      display: grid;
      grid-template-columns: 120px 1fr;
      gap: 8px 12px;
      margin: 0;
      font-size: 15px;
    }
    #mackey-visualizer .mv-panel dt {
      font-weight: 700;
    }
    #mackey-visualizer .mv-panel dd {
      margin: 0;
      overflow-wrap: anywhere;
    }
    #mackey-visualizer code {
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 13px;
      background: #f6f6f6;
      border: 1px solid #000;
      padding: 2px 5px;
      white-space: pre-wrap;
      word-break: break-word;
    }
    #mackey-visualizer .mv-code-line {
      display: flex;
      flex-wrap: wrap;
      gap: 7px;
      align-items: baseline;
      margin: 8px 0 0;
    }
    #mackey-visualizer .mv-code-label {
      font-size: 14px;
      font-weight: 700;
    }
    #mackey-visualizer .mv-code-value {
      font-size: 13px;
    }
    #mackey-visualizer .mv-math {
      margin: 8px 0 14px;
      padding: 10px 12px;
      border: 1px solid #000;
      background: #fff;
      overflow-x: auto;
    }
    #mackey-visualizer .mv-math .katex-display {
      margin: 0;
    }
    #mackey-visualizer .mv-math-inline {
      margin: 0;
      padding: 0;
      border: none;
      display: inline-block;
      overflow: visible;
    }
    #mackey-visualizer .mv-summary {
      margin: 8px 0;
      font-size: 15px;
    }
    #mackey-visualizer .node rect {
      fill: #fff;
      stroke: #000;
      stroke-width: 1.8;
      rx: 6;
      ry: 6;
    }
    #mackey-visualizer .node:hover rect,
    #mackey-visualizer .node.selected rect {
      stroke-width: 3;
    }
    #mackey-visualizer .mv-node-box {
      width: 172px;
      height: 78px;
      display: flex;
      flex-direction: column;
      align-items: center;
      justify-content: center;
      gap: 4px;
      text-align: center;
      color: #000;
      font-size: 16px;
      line-height: 1.15;
      overflow: hidden;
      pointer-events: none;
    }
    #mackey-visualizer .mv-node-box .katex {
      font-size: 1em;
    }
    #mackey-visualizer .edge {
      fill: none;
      stroke: #000;
      stroke-width: 2.6;
      stroke-linecap: round;
      cursor: pointer;
      pointer-events: none;
    }
    #mackey-visualizer .edge.selected,
    #mackey-visualizer .conj.selected {
      stroke-width: 4.4;
    }
    #mackey-visualizer .edge-hit {
      fill: none;
      stroke: transparent;
      stroke-width: 20;
      cursor: pointer;
      pointer-events: stroke;
    }
    #mackey-visualizer .spine {
      fill: none;
      stroke: #000;
      stroke-width: 1.4;
      stroke-dasharray: 6 6;
      opacity: 0.22;
      pointer-events: none;
    }
    #mackey-visualizer .conj {
      fill: none;
      stroke: #000;
      stroke-width: 2.4;
      stroke-dasharray: 4 4;
      stroke-linecap: round;
      pointer-events: none;
    }
  </style>
  <div class="mv-toolbar">
    <div class="mv-toolbar-group">
      <button data-kind="res" class="active">Res</button>
      <button data-kind="tr" class="active">Tr</button>
      <button data-kind="conj">Conj</button>
    </div>
    <div class="mv-toolbar-group">
      <button data-action="fit" class="active">Fit width</button>
      <button data-action="expand">Expand graph</button>
    </div>
  </div>
  <div class="mv-shell">
    <div class="mv-graph">
      <svg viewBox="0 0 900 620" role="img" aria-label="Mackey functor subgroup lattice"></svg>
    </div>
    <aside class="mv-panel open">
      <h3>Selection</h3>
      <div class="panel-content">Click a subgroup node.</div>
    </aside>
  </div>
</div>
<script>
(function() {
  const root = document.currentScript.previousElementSibling;
  const data = $data_json;
  const graph = root.querySelector(".mv-graph");
  const svg = root.querySelector("svg");
  const panelBox = root.querySelector(".mv-panel");
  const panel = root.querySelector(".panel-content");
  const buttons = root.querySelectorAll("button[data-kind]");
  const actionButtons = root.querySelectorAll("button[data-action]");
  const state = { res: true, tr: true, conj: false, selected: [], selectedEdge: null, fit: true, expanded: false };
  const baseWidth = 900;
  const margin = { left: 104, right: 104, top: 92, bottom: 92 };
  const nodeWidth = 172;
  const nodeHeight = 78;
  let viewWidth = baseWidth;
  let viewHeight = 620;
  let katexRequested = false;

  function esc(x) {
    return String(x).replace(/[&<>"']/g, function(c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }

  function ensureKatex() {
    if (window.katex || katexRequested) return;
    katexRequested = true;
    if (!document.getElementById("mackey-katex-css")) {
      const link = document.createElement("link");
      link.id = "mackey-katex-css";
      link.rel = "stylesheet";
      link.href = "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css";
      document.head.appendChild(link);
    }
    const script = document.createElement("script");
    script.src = "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js";
    script.async = true;
    script.onload = function() { draw(); };
    document.head.appendChild(script);
  }

  function math(tex, displayMode) {
    const isDisplay = displayMode === true;
    if (window.katex) {
      return window.katex.renderToString(String(tex), { displayMode: isDisplay, throwOnError: false });
    }
    ensureKatex();
    return "<span class='tex-source'>\\\\(" + esc(tex) + "\\\\)</span>";
  }

  function code(x) {
    return "<code>" + esc(x) + "</code>";
  }

  function codeLine(label, value) {
    if (value === null || value === undefined) return "";
    const text = String(value).trim();
    if (!text) return "";
    return "<div class='mv-code-line'><span class='mv-code-label'>" + esc(label) + "</span><code class='mv-code-value'>" + esc(text) + "</code></div>";
  }

  function matrixTex(matrix) {
    if (!matrix || matrix.length === 0) return "0";
    return "\\\\begin{bmatrix}" + matrix.map(function(row) {
      return row.map(function(x) { return String(x); }).join(" & ");
    }).join(" \\\\\\\\ ") + "\\\\end{bmatrix}";
  }

  function matrixHTML(matrix) {
    return "<div class='mv-math'>" + math(matrixTex(matrix), true) + "</div>";
  }

  function byNode(id) {
    return data.nodes.find(function(n) { return n.id === id; });
  }

  function matrixBetween(entries, source, target) {
    return entries.find(function(e) { return e.source === source && e.target === target; });
  }

  function comparable(a, b) {
    return matrixBetween(data.restrictions, a, b) ||
           matrixBetween(data.restrictions, b, a) ||
           matrixBetween(data.transfers, a, b) ||
           matrixBetween(data.transfers, b, a);
  }

  function layoutNodes() {
    const byId = new Map(data.nodes.map(function(n) { return [n.id, n]; }));
    const children = new Map(data.nodes.map(function(n) { return [n.id, []]; }));
    const indegree = new Map(data.nodes.map(function(n) { return [n.id, 0]; }));

    data.covers.forEach(function(e) {
      const large = e.target;
      const small = e.source;
      children.get(large).push(small);
      indegree.set(small, indegree.get(small) + 1);
    });

    const rank = new Map();
    const queue = data.nodes
      .filter(function(n) { return indegree.get(n.id) === 0; })
      .sort(function(a, b) { return b.size - a.size || a.id - b.id; })
      .map(function(n) { rank.set(n.id, 0); return n.id; });

    while (queue.length > 0) {
      const id = queue.shift();
      const nextRank = rank.get(id) + 1;
      children.get(id).forEach(function(child) {
        rank.set(child, Math.max(rank.get(child) || 0, nextRank));
        indegree.set(child, indegree.get(child) - 1);
        if (indegree.get(child) === 0) queue.push(child);
      });
    }

    data.nodes.forEach(function(n) {
      if (!rank.has(n.id)) rank.set(n.id, 0);
    });

    const levels = {};
    data.nodes.forEach(function(n) {
      const key = String(rank.get(n.id));
      if (!levels[key]) levels[key] = [];
      levels[key].push(n);
    });

    const ranks = Object.keys(levels).map(Number).sort(function(a, b) { return a - b; });
    const maxLevelSize = Math.max.apply(null, ranks.map(function(r) { return levels[String(r)].length; }));
    const levelGap = Math.min(320, Math.max(210, 178 + Math.max(maxLevelSize - 1, 0) * 20));
    viewWidth = Math.max(baseWidth, margin.left + margin.right + maxLevelSize * nodeWidth + Math.max(maxLevelSize - 1, 0) * 104);
    viewHeight = margin.top + margin.bottom + Math.max(ranks.length - 1, 1) * levelGap;
    svg.setAttribute("viewBox", "0 0 " + viewWidth + " " + viewHeight);
    svg.style.width = state.fit ? "100%" : viewWidth + "px";
    svg.style.height = viewHeight + "px";

    ranks.forEach(function(r, levelIndex) {
      const nodes = levels[String(r)].sort(function(a, b) {
        return a.label.localeCompare(b.label) || a.id - b.id;
      });
      const y = ranks.length === 1
        ? viewHeight / 2
        : margin.top + levelIndex * ((viewHeight - margin.top - margin.bottom) / (ranks.length - 1));
      nodes.forEach(function(n, i) {
        const x = nodes.length === 1
          ? viewWidth / 2
          : margin.left + i * ((viewWidth - margin.left - margin.right) / (nodes.length - 1));
        n.x = x;
        n.y = y;
      });
    });
  }

  function arrowMarker(id, dashed) {
    return "<marker id='" + id + "' markerWidth='10' markerHeight='10' refX='9' refY='5' orient='auto' markerUnits='strokeWidth'>" +
      "<path d='M0,0 L10,5 L0,10 Z' fill='#000'></path></marker>";
  }

  function edgePoint(node, toward) {
    const dx = toward.x - node.x;
    const dy = toward.y - node.y;
    const scale = Math.max(Math.abs(dx) / (nodeWidth / 2), Math.abs(dy) / (nodeHeight / 2), 1);
    return { x: node.x + dx / scale, y: node.y + dy / scale };
  }

  function curve(aNode, bNode, offset) {
    const a = edgePoint(aNode, bNode);
    const b = edgePoint(bNode, aNode);
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const len = Math.max(Math.sqrt(dx * dx + dy * dy), 1);
    const nx = -dy / len;
    const ny = dx / len;
    const mx = (a.x + b.x) / 2 + nx * offset;
    const my = (a.y + b.y) / 2 + ny * offset;
    return "M" + a.x + "," + a.y + " Q" + mx + "," + my + " " + b.x + "," + b.y;
  }

  function selfLoop(n) {
    const x = n.x + nodeWidth / 2 - 12;
    const y = n.y - nodeHeight / 2 + 4;
    return "M" + x + "," + y + " C" + (x + 46) + "," + (y - 50) + " " + (x + 82) + "," + (y + 6) + " " + (x + 20) + "," + (y + 20);
  }

  function draw() {
    layoutNodes();
    let html = "<defs>" + arrowMarker("arrow", false) + arrowMarker("arrow-dashed", true) + "</defs>";

    data.covers.forEach(function(e) {
      const small = byNode(e.source);
      const large = byNode(e.target);
      html += "<path class='spine' d='" + curve(large, small, 0) + "'></path>";
    });

    if (state.res) {
      data.covers.forEach(function(e) {
        const small = byNode(e.source);
        const large = byNode(e.target);
        const idx = data.restrictions.findIndex(function(m) { return m.source === large.id && m.target === small.id; });
        const selected = state.selectedEdge && state.selectedEdge.kind === "res" && state.selectedEdge.index === idx ? " selected" : "";
        const path = curve(large, small, -30);
        html += "<path class='edge-hit clickedge' data-kind='res' data-index='" + idx + "' d='" + path + "'></path>";
        html += "<path class='edge res" + selected + "' marker-end='url(#arrow)' d='" + path + "'></path>";
      });
    }

    if (state.tr) {
      data.covers.forEach(function(e) {
        const small = byNode(e.source);
        const large = byNode(e.target);
        const idx = data.transfers.findIndex(function(m) { return m.source === small.id && m.target === large.id; });
        const selected = state.selectedEdge && state.selectedEdge.kind === "tr" && state.selectedEdge.index === idx ? " selected" : "";
        const path = curve(small, large, -30);
        html += "<path class='edge-hit clickedge' data-kind='tr' data-index='" + idx + "' d='" + path + "'></path>";
        html += "<path class='edge tr" + selected + "' marker-end='url(#arrow)' d='" + path + "'></path>";
      });
    }

    if (state.conj) {
      data.conjugations.forEach(function(e, idx) {
        const source = byNode(e.source);
        const target = byNode(e.target);
        const path = source.id === target.id ? selfLoop(source) : curve(source, target, 32);
        const selected = state.selectedEdge && state.selectedEdge.kind === "conj" && state.selectedEdge.index === idx ? " selected" : "";
        html += "<path class='edge-hit clickedge' data-kind='conj' data-index='" + idx + "' d='" + path + "'></path>";
        html += "<path class='conj" + selected + "' marker-end='url(#arrow-dashed)' d='" + path + "'></path>";
      });
    }

    data.nodes.forEach(function(n) {
      const selected = state.selected.indexOf(n.id) !== -1 ? " selected" : "";
      html += "<g class='node" + selected + "' data-id='" + n.id + "' transform='translate(" + (n.x - nodeWidth / 2) + "," + (n.y - nodeHeight / 2) + ")'>";
      html += "<rect width='" + nodeWidth + "' height='" + nodeHeight + "' rx='6' ry='6'></rect>";
      html += "<foreignObject width='" + nodeWidth + "' height='" + nodeHeight + "'>";
      html += "<div xmlns='http://www.w3.org/1999/xhtml' class='mv-node-box'>";
      html += "<div>" + math(n.label_tex || n.label) + "</div>";
      html += "<div>" + math(n.value_tex || n.value) + "</div>";
      html += "</div></foreignObject>";
      html += "</g>";
    });

    svg.innerHTML = html;
    svg.querySelectorAll(".node").forEach(function(node) {
      node.addEventListener("click", function(event) {
        event.stopPropagation();
        selectNode(Number(node.getAttribute("data-id")));
      });
    });
    svg.querySelectorAll(".clickedge").forEach(function(edge) {
      edge.addEventListener("click", function(event) {
        event.stopPropagation();
        selectEdge(edge.getAttribute("data-kind"), Number(edge.getAttribute("data-index")));
      });
    });
  }

  function nodeDetails(id) {
    const n = byNode(id);
    const conjs = data.conjugations.filter(function(e) { return e.source === id; });
    let html = "<h3>Subgroup</h3>";
    html += "<div class='mv-math'>" + math(n.orbit_tex || n.orbit) + "</div>";
    html += "<div class='mv-summary'><strong>Value</strong> <span class='mv-math-inline'>" + math(n.value_tex || n.value) + "</span></div>";
    html += "<div class='mv-summary'><strong>Subgroup order</strong> <span class='mv-math-inline'>" + math(String(n.size)) + "</span></div>";
    html += codeLine("Lattice id", n.id);
    html += codeLine("GAP subgroup object", n.gap_name);
    html += codeLine("GAP subgroup IdGroup", n.id_group);
    html += codeLine("GAP subgroup generators", n.generators);
    html += "<h4>Generator conjugations</h4>";
    if (conjs.length === 0) {
      html += "<p>None stored for group generators.</p>";
    } else {
      conjs.forEach(function(c) {
        html += "<div><strong>" + math((n.label_tex || n.label) + " \\\\xrightarrow{" + esc(c.element) + "} " + (byNode(c.target).label_tex || byNode(c.target).label)) + "</strong>";
        html += matrixHTML(c.matrix) + "</div>";
      });
    }
    openPanel(html);
  }

  function pairDetails(a, b) {
    const A = byNode(a);
    const B = byNode(b);
    const resAB = matrixBetween(data.restrictions, a, b);
    const resBA = matrixBetween(data.restrictions, b, a);
    const trAB = matrixBetween(data.transfers, a, b);
    const trBA = matrixBetween(data.transfers, b, a);

    let html = "<h3>Comparable Pair</h3><div class='mv-summary'><strong>" + math(A.label_tex || A.label) + "</strong> and <strong>" + math(B.label_tex || B.label) + "</strong></div>";

    if (!comparable(a, b)) {
      openPanel(html + "<p>These subgroups are not comparable in the lattice.</p>");
      return;
    }

    if (resAB) html += "<h4>Res: " + math((A.label_tex || A.label) + " \\\\to " + (B.label_tex || B.label)) + "</h4>" + matrixHTML(resAB.matrix);
    if (resBA) html += "<h4>Res: " + math((B.label_tex || B.label) + " \\\\to " + (A.label_tex || A.label)) + "</h4>" + matrixHTML(resBA.matrix);
    if (trAB) html += "<h4>Tr: " + math((A.label_tex || A.label) + " \\\\to " + (B.label_tex || B.label)) + "</h4>" + matrixHTML(trAB.matrix);
    if (trBA) html += "<h4>Tr: " + math((B.label_tex || B.label) + " \\\\to " + (A.label_tex || A.label)) + "</h4>" + matrixHTML(trBA.matrix);
    openPanel(html);
  }

  function edgeDetails(kind, entry) {
    const source = byNode(entry.source);
    const target = byNode(entry.target);
    const title = kind === "res" ? "Restriction" : kind === "tr" ? "Transfer" : "Conjugation";
    let html = "<h3>" + title + "</h3>";
    html += "<div class='mv-math'>" + math((source.label_tex || source.label) + " \\\\to " + (target.label_tex || target.label)) + "</div>";
    if (kind === "conj") {
      html += codeLine("Conjugating element", entry.element);
    }
    html += matrixHTML(entry.matrix);
    openPanel(html);
  }

  function openPanel(html) {
    panel.innerHTML = html;
    panelBox.classList.add("open");
  }

  function clearSelection() {
    state.selected = [];
    state.selectedEdge = null;
    panel.innerHTML = "Click a subgroup node or map arrow.";
    panelBox.classList.remove("open");
  }

  function selectEdge(kind, index) {
    if (index < 0) return;
    state.selected = [];
    state.selectedEdge = { kind: kind, index: index };
    const entry = kind === "res"
      ? data.restrictions[index]
      : kind === "tr"
        ? data.transfers[index]
        : data.conjugations[index];
    edgeDetails(kind, entry);
    draw();
  }

  function selectNode(id) {
    state.selectedEdge = null;
    if (state.selected.length === 0 || state.selected[0] === id) {
      state.selected = [id];
      nodeDetails(id);
    } else {
      state.selected = [state.selected[0], id];
      pairDetails(state.selected[0], id);
    }
    draw();
  }

  buttons.forEach(function(button) {
    button.addEventListener("click", function() {
      const kind = button.getAttribute("data-kind");
      state[kind] = !state[kind];
      button.classList.toggle("active", state[kind]);
      draw();
    });
  });

  actionButtons.forEach(function(button) {
    button.addEventListener("click", function() {
      const action = button.getAttribute("data-action");
      if (action === "fit") {
        state.fit = !state.fit;
        button.classList.toggle("active", state.fit);
        draw();
      } else if (action === "expand") {
        state.expanded = !state.expanded;
        root.classList.toggle("expanded", state.expanded);
        button.classList.toggle("active", state.expanded);
        button.textContent = state.expanded ? "Collapse graph" : "Expand graph";
        draw();
      }
    });
  });

  svg.addEventListener("click", function() {
    clearSelection();
    draw();
  });

  draw();
})();
</script>
""")
end

# ╔═╡ 55555555-5555-5555-5555-555555555555
mackey_visualizer_html(M)

# ╔═╡ Cell order:
# ╠═11111111-1111-1111-1111-111111111111
# ╟─22222222-2222-2222-2222-222222222222
# ╠═33333333-3333-3333-3333-333333333333
# ╠═44444444-4444-4444-4444-444444444444
# ╠═55555555-5555-5555-5555-555555555555

begin
    html = repr(MIME"text/html"(), mackey_visualizer_html(M))
    output_path = joinpath(@__DIR__, "mackey_functor_visualizer.html")

    write(output_path, """
    <!doctype html>
    <html>
    <head>
      <meta charset="utf-8">
      <title>Mackey Functor Visualizer</title>
    </head>
    <body>
    $html
    </body>
    </html>
    """)

    output_path
end