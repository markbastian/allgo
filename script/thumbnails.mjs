// The gallery's thumbnails: a screenshot of every demo in demos.json,
// taken from a running `make serve` by a headless Chrome.
//
//   make thumbs                      (with `make serve` running)
//   BASE=http://localhost:3000/ CHROME=/path/to/chrome node script/thumbnails.mjs [id ...]
//
// Each demo is opened with `?thumb`, which hides its title bar, controls
// and readouts (see allgo.app), left to run for a few seconds of real time
// so there is something to see -- a flock that has flocked, a fire that
// has caught -- and captured at 1280x720, scaled to 640x360 WebP, into
// resources/public/thumbs/<id>.webp. Real time rather than headless
// Chrome's virtual time: under that the simulations barely advance, since
// a frame is only as fast as the software renderer can draw it.
//
// Chrome is driven over the DevTools protocol with the `ws` package that
// is already in node_modules; nothing else is needed.

import { spawn } from "node:child_process";
import { mkdtempSync, readFileSync, writeFileSync, mkdirSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { createRequire } from "node:module";

const require = createRequire(import.meta.url);
const WebSocket = require("ws");

const BASE = process.env.BASE || "http://localhost:3000/";
const CHROME = process.env.CHROME ||
  "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";
const PORT = 9333;
const OUT = "resources/public/thumbs";

// How long each demo runs before its picture is taken, where the default
// is too soon: the ones that build something up over time.
const SETTLE = { fire: 12000, dla: 15000, erosion: 12000, flip: 9000,
                 "sphere-fluid": 10000, island: 9000, planet: 12000,
                 "boids-voronoi-3d": 8000, "orbit-determination": 8000 };
const DEFAULT_SETTLE = 6000;

// How far to zoom toward the middle, for the demos whose subject is small
// in a wide window: their cameras were framed for the square the old demo
// page drew them in. 2 keeps the middle half, and so on.
const ZOOM = { "soft-body": 2.2, skinning: 2.2, cloth: 1.7, joints: 1.8,
               rigid: 1.6, "human-arm": 1.4, kepler: 1.8, hex: 1.4,
               bricks: 1.9, ragdoll: 1.7 };

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

const catalog = JSON.parse(readFileSync("resources/public/demos.json", "utf8"));
const only = process.argv.slice(2);
const demos = catalog.demos.filter((d) => only.length === 0 || only.includes(d.id));

const profile = mkdtempSync(join(tmpdir(), "allgo-thumbs-"));
const chrome = spawn(CHROME, [
  "--headless=new", `--remote-debugging-port=${PORT}`, `--user-data-dir=${profile}`,
  "--window-size=1280,720", "--hide-scrollbars", "--mute-audio",
  // Software WebGL, so this works on a machine with no GPU to hand.
  "--use-angle=swiftshader", "--enable-unsafe-swiftshader",
  "about:blank",
], { stdio: "ignore" });

let exited = null;
chrome.on("exit", (code) => { exited = code; });

async function target() {
  let last = null;
  for (let i = 0; i < 150; i++) {
    if (exited !== null) throw new Error(`Chrome exited (${exited}) before it could be reached`);
    try {
      const list = await (await fetch(`http://127.0.0.1:${PORT}/json/list`)).json();
      // The tab, not an extension's background page or the omnibox.
      const page = list.find((t) => t.type === "page");
      if (page) return page.webSocketDebuggerUrl;
    } catch (e) { last = e; }
    await sleep(200);
  }
  throw new Error(`Chrome did not answer on port ${PORT}: ${last && (last.cause || last).message}`);
}

const ws = new WebSocket(await target(), { perMessageDeflate: false });
await new Promise((r) => ws.once("open", r));
let nextId = 0;
const pending = new Map();
const waiters = [];
ws.on("message", (raw) => {
  const msg = JSON.parse(raw);
  if (msg.id !== undefined && pending.has(msg.id)) {
    const { resolve, reject } = pending.get(msg.id);
    pending.delete(msg.id);
    msg.error ? reject(new Error(msg.error.message)) : resolve(msg.result);
  } else if (msg.method) {
    waiters.filter((w) => w.method === msg.method).forEach((w) => w.resolve(msg.params));
  }
});
const send = (method, params = {}) => new Promise((resolve, reject) => {
  const id = ++nextId;
  pending.set(id, { resolve, reject });
  ws.send(JSON.stringify({ id, method, params }));
});
const once = (method, ms) => new Promise((resolve) => {
  const w = { method, resolve };
  waiters.push(w);
  setTimeout(resolve, ms);
});

await send("Page.enable");
await send("Emulation.setDeviceMetricsOverride",
           { width: 1280, height: 720, deviceScaleFactor: 1, mobile: false });
mkdirSync(OUT, { recursive: true });

for (const d of demos) {
  const url = new URL(d.href, BASE);
  url.searchParams.set("thumb", "");
  const loaded = once("Page.loadEventFired", 20000);
  await send("Page.navigate", { url: url.href });
  await loaded;
  await sleep(SETTLE[d.id] || DEFAULT_SETTLE);
  const z = ZOOM[d.id] || 1;
  const w = 1280 / z, h = 720 / z;
  const { data } = await send("Page.captureScreenshot", {
    format: "webp", quality: 78,
    clip: { x: (1280 - w) / 2, y: (720 - h) / 2, width: w, height: h, scale: 640 / w },
  });
  const file = join(OUT, `${d.id}.webp`);
  writeFileSync(file, Buffer.from(data, "base64"));
  console.log(`${d.id.padEnd(22)} ${(data.length * 0.75 / 1024).toFixed(0)} KB`);
}

ws.close();
chrome.kill();
rmSync(profile, { recursive: true, force: true });
