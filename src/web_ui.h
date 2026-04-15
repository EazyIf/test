#ifndef WEB_UI_H
#define WEB_UI_H

#include <string>

static const std::string HTML_PAGE = R"HTML(<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>3D Aerodynamics Simulator</title>
<style>
* { margin:0; padding:0; box-sizing:border-box; }
body { background:#1a1a2e; color:#eee; font-family:'Segoe UI',system-ui,sans-serif; overflow:hidden; }
#canvas3d { position:fixed; top:0; left:0; width:100%; height:100%; }

#hud {
  position:fixed; top:16px; left:16px; z-index:10;
  background:rgba(10,10,30,0.85); border:1px solid #334; border-radius:12px;
  padding:16px 20px; min-width:280px; backdrop-filter:blur(8px);
}
#hud h1 { font-size:15px; color:#6cf; margin-bottom:8px; font-weight:600; }
#hud .stat { font-size:13px; margin:3px 0; color:#bbb; }
#hud .stat b { color:#fff; font-variant-numeric:tabular-nums; }
#hud .force-row { display:flex; gap:18px; margin:6px 0; }
#hud .force-val { font-size:22px; font-weight:700; }
#hud .force-val.cd { color:#ff6b6b; }
#hud .force-val.cl { color:#51cf66; }
#hud .force-val.cs { color:#74c0fc; }
#hud .force-label { font-size:10px; text-transform:uppercase; color:#888; letter-spacing:1px; }

#controls {
  position:fixed; top:16px; right:16px; z-index:10;
  background:rgba(10,10,30,0.85); border:1px solid #334; border-radius:12px;
  padding:16px 20px; width:260px; backdrop-filter:blur(8px);
}
#controls h2 { font-size:13px; color:#6cf; margin-bottom:10px; }
.ctrl-row { margin:8px 0; }
.ctrl-row label { font-size:11px; color:#999; display:block; margin-bottom:2px; }
.ctrl-row input[type=range] { width:100%; accent-color:#6cf; }
.ctrl-row .val { float:right; color:#fff; font-size:11px; font-variant-numeric:tabular-nums; }
.btn-row { display:flex; gap:8px; margin-top:12px; }
.btn { background:#334; color:#ccc; border:1px solid #556; border-radius:6px;
       padding:6px 12px; font-size:12px; cursor:pointer; flex:1; text-align:center; }
.btn:hover { background:#445; color:#fff; }
.btn.active { background:#264; border-color:#4a8; color:#6f6; }

#status { position:fixed; bottom:16px; left:50%; transform:translateX(-50%); z-index:10;
  background:rgba(10,10,30,0.8); border:1px solid #334; border-radius:8px;
  padding:8px 16px; font-size:11px; color:#888; }
#legend { position:fixed; bottom:16px; right:16px; z-index:10;
  background:rgba(10,10,30,0.8); border:1px solid #334; border-radius:8px;
  padding:8px 12px; font-size:11px; color:#888; }
</style>
</head>
<body>
<canvas id="canvas3d"></canvas>

<div id="hud">
  <h1>3D Aerodynamics Simulator</h1>
  <div class="stat">Step: <b id="s_step">0</b> &nbsp; Vehicle: <b id="s_vehicle">-</b></div>
  <div class="stat">Grid: <b id="s_grid">-</b> &nbsp; tau: <b id="s_tau">-</b> &nbsp; Ma: <b id="s_mach">-</b></div>
  <hr style="border-color:#333;margin:8px 0">
  <div class="force-row">
    <div><div class="force-label">Drag (Cd)</div><div class="force-val cd" id="s_cd">-</div></div>
    <div><div class="force-label">Lift (Cl)</div><div class="force-val cl" id="s_cl">-</div></div>
    <div><div class="force-label">Side (Cs)</div><div class="force-val cs" id="s_cs">-</div></div>
  </div>
  <div class="stat">Residual: <b id="s_res">-</b> &nbsp; V<sub>max</sub>: <b id="s_vmax">-</b></div>
</div>

<div id="controls">
  <h2>Parameters</h2>
  <div class="ctrl-row">
    <label>Reynolds <span class="val" id="v_re">150</span></label>
    <input type="range" id="r_re" min="25" max="500" step="25" value="150">
  </div>
  <div class="ctrl-row" id="row_aoa" style="display:none">
    <label>Angle of Attack <span class="val" id="v_aoa">5</span>&deg;</label>
    <input type="range" id="r_aoa" min="-10" max="20" step="1" value="5">
  </div>
  <div class="ctrl-row">
    <label>Yaw <span class="val" id="v_yaw">0</span>&deg;</label>
    <input type="range" id="r_yaw" min="-30" max="30" step="1" value="0">
  </div>
  <div class="ctrl-row" id="row_slant">
    <label>Slant Angle <span class="val" id="v_slant">25</span>&deg;</label>
    <input type="range" id="r_slant" min="5" max="40" step="1" value="25">
  </div>
  <div class="ctrl-row" id="row_sail" style="display:none">
    <label>Sail Angle <span class="val" id="v_sail">15</span>&deg;</label>
    <input type="range" id="r_sail" min="-45" max="45" step="1" value="15">
  </div>
  <div class="btn-row">
    <div class="btn" id="btn_pause" onclick="togglePause()">Pause</div>
    <div class="btn" onclick="sendCmd('reset')">Reset</div>
  </div>
  <div class="btn-row">
    <div class="btn" onclick="sendCmd('vehicle')">Switch Vehicle</div>
  </div>
</div>

<div id="status">Connecting to simulation...</div>
<div id="legend">Drag to rotate &bull; Scroll to zoom &bull; Right-drag to pan</div>

<script type="importmap">
{ "imports": {
    "three": "https://cdn.jsdelivr.net/npm/three@0.164.1/build/three.module.js",
    "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.164.1/examples/jsm/"
}}
</script>
<script type="module">
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

// --- Scene setup ---
const canvas = document.getElementById('canvas3d');
const renderer = new THREE.WebGLRenderer({ canvas, antialias:true });
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
renderer.setClearColor(0x0a0a1e);
renderer.toneMapping = THREE.ACESFilmicToneMapping;

const scene = new THREE.Scene();
scene.fog = new THREE.FogExp2(0x0a0a1e, 0.003);

const camera = new THREE.PerspectiveCamera(55, window.innerWidth/window.innerHeight, 0.1, 2000);
camera.position.set(80, 50, 60);

const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
controls.dampingFactor = 0.08;
controls.target.set(40, 0, 0);

// Lights
scene.add(new THREE.AmbientLight(0x334466, 1.5));
const dirLight = new THREE.DirectionalLight(0xffffff, 2.0);
dirLight.position.set(50, 80, 40);
scene.add(dirLight);
const dirLight2 = new THREE.DirectionalLight(0x4488cc, 0.8);
dirLight2.position.set(-30, -20, 60);
scene.add(dirLight2);

// Grid floor
const gridHelper = new THREE.GridHelper(200, 40, 0x222244, 0x111133);
gridHelper.rotation.x = Math.PI / 2;
scene.add(gridHelper);

// Wind direction arrow
const windArrow = new THREE.ArrowHelper(
  new THREE.Vector3(1,0,0), new THREE.Vector3(-30,0,0), 20, 0x4488ff, 4, 2);
scene.add(windArrow);

// --- Body mesh ---
let bodyMesh = null;
let arrowGroup = new THREE.Group();
scene.add(arrowGroup);

function buildBodyMesh(data) {
  if (bodyMesh) { scene.remove(bodyMesh); bodyMesh.geometry.dispose(); bodyMesh.material.dispose(); }

  const cells = data.body;
  if (!cells || cells.length === 0) return;

  const count = cells.length;
  const geo = new THREE.BoxGeometry(1.0, 1.0, 1.0);
  const mat = new THREE.MeshPhysicalMaterial({
    vertexColors: true, roughness: 0.4, metalness: 0.1,
    clearcoat: 0.3, clearcoatRoughness: 0.2
  });

  const instGeo = new THREE.InstancedBufferGeometry().copy(geo);
  // Use InstancedMesh
  bodyMesh = new THREE.InstancedMesh(geo, mat, count);

  const dummy = new THREE.Object3D();
  const color = new THREE.Color();

  // Pressure range for coloring
  let pMin = Infinity, pMax = -Infinity;
  for (const c of cells) {
    if (c[6] < pMin) pMin = c[6];
    if (c[6] > pMax) pMax = c[6];
  }
  const pRange = Math.max(pMax - pMin, 1e-10);

  const cx = data.cx || 0, cy = data.cy || 0, cz = data.cz || 0;

  for (let i = 0; i < count; i++) {
    const [x, y, z, nx, ny, nz, pressure] = cells[i];
    dummy.position.set(x - cx, y - cy, z - cz);
    dummy.updateMatrix();
    bodyMesh.setMatrixAt(i, dummy.matrix);

    // Cool-to-hot pressure coloring
    const t = (pressure - pMin) / pRange;
    color.setHSL(0.65 - t * 0.65, 0.85, 0.45 + t * 0.2);
    bodyMesh.setColorAt(i, color);
  }

  bodyMesh.instanceMatrix.needsUpdate = true;
  bodyMesh.instanceColor.needsUpdate = true;
  scene.add(bodyMesh);

  // Center camera on body
  controls.target.set(0, 0, 0);
}

function buildFlowArrows(data) {
  // Remove old arrows
  while (arrowGroup.children.length > 0) {
    const c = arrowGroup.children[0];
    arrowGroup.remove(c);
    if (c.geometry) c.geometry.dispose();
    if (c.material) c.material.dispose();
  }

  const flow = data.flow;
  if (!flow || flow.length === 0) return;

  const cx = data.cx || 0, cy = data.cy || 0, cz = data.cz || 0;

  // Use line segments for performance
  const positions = [];
  const colors = [];
  const color = new THREE.Color();

  const vMax = data.max_vel || 0.1;

  for (const f of flow) {
    const [x, y, z, vx, vy, vz] = f;
    const vmag = Math.sqrt(vx*vx + vy*vy + vz*vz);
    if (vmag < 1e-8) continue;

    const px = x - cx, py = y - cy, pz = z - cz;
    const scale = 15.0 / Math.max(vMax, 0.01);

    positions.push(px, py, pz);
    positions.push(px + vx*scale, py + vy*scale, pz + vz*scale);

    const t = Math.min(vmag / Math.max(vMax, 0.01), 1.0);
    color.setHSL(0.6 - t * 0.6, 0.9, 0.4 + t * 0.3);
    colors.push(color.r, color.g, color.b);
    colors.push(color.r, color.g, color.b);
  }

  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  geo.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
  const mat = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.6 });
  arrowGroup.add(new THREE.LineSegments(geo, mat));
}

// --- State polling ---
let paused = false;
let meshLoaded = false;
let lastMeshVehicle = '';

async function fetchState() {
  try {
    const res = await fetch('/state');
    const d = await res.json();

    document.getElementById('s_step').textContent = d.step;
    document.getElementById('s_vehicle').textContent = d.vehicle;
    document.getElementById('s_grid').textContent = d.nx+'x'+d.ny+'x'+d.nz;
    document.getElementById('s_tau').textContent = d.tau.toFixed(4);
    document.getElementById('s_mach').textContent = (d.max_vel/0.5774).toFixed(3);
    document.getElementById('s_cd').textContent = d.cd.toFixed(4);
    document.getElementById('s_cl').textContent = d.cl.toFixed(4);
    document.getElementById('s_cs').textContent = d.cs.toFixed(4);
    document.getElementById('s_res').textContent = d.residual.toExponential(2);
    document.getElementById('s_vmax').textContent = d.max_vel.toFixed(4);

    document.getElementById('status').textContent =
      d.paused ? 'PAUSED' : d.unstable ? 'UNSTABLE - lower Re or reset' : 'Step ' + d.step;

    // Show/hide vehicle-specific controls
    const v = d.vehicle;
    document.getElementById('row_aoa').style.display = (v === 'plane') ? '' : 'none';
    document.getElementById('row_slant').style.display = (v === 'car') ? '' : 'none';
    document.getElementById('row_sail').style.display = (v === 'sailboat') ? '' : 'none';

    // Reload mesh if vehicle changed
    if (v !== lastMeshVehicle || !meshLoaded) {
      lastMeshVehicle = v;
      await fetchMesh();
    }

    // Update flow arrows
    buildFlowArrows(d);

  } catch(e) {
    document.getElementById('status').textContent = 'Connection lost - retrying...';
  }
}

async function fetchMesh() {
  try {
    const res = await fetch('/mesh');
    const d = await res.json();
    buildBodyMesh(d);
    meshLoaded = true;
  } catch(e) {}
}

// --- Controls ---
function sendParam(name, value) {
  fetch('/control', {
    method: 'POST',
    body: JSON.stringify({action:'param', name, value})
  });
}

function sendCmd(cmd) {
  fetch('/control', { method:'POST', body: JSON.stringify({action:cmd}) });
  if (cmd === 'vehicle') { meshLoaded = false; }
}

window.sendCmd = sendCmd;

function togglePause() {
  paused = !paused;
  sendCmd('pause');
  document.getElementById('btn_pause').textContent = paused ? 'Resume' : 'Pause';
  document.getElementById('btn_pause').classList.toggle('active', paused);
}
window.togglePause = togglePause;

// Wire up sliders
for (const [id, name] of [['r_re','re'],['r_aoa','aoa'],['r_yaw','yaw'],['r_slant','slant'],['r_sail','sail']]) {
  const el = document.getElementById(id);
  const vEl = document.getElementById('v_'+name);
  el.addEventListener('input', () => {
    vEl.textContent = el.value;
    sendParam(name, parseFloat(el.value));
  });
}

// --- Render loop ---
function animate() {
  requestAnimationFrame(animate);
  controls.update();
  renderer.render(scene, camera);
}
animate();

// Poll state
setInterval(fetchState, 250);
setTimeout(fetchMesh, 500);

window.addEventListener('resize', () => {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
});
</script>
</body>
</html>
)HTML";

#endif
