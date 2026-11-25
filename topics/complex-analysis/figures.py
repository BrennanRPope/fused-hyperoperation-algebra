"""Domain colouring of the lift 1 + z: hue = angle seen from the wall, bands = log distance.
Writes figures/phase.png, phase.jpg and figures/phase_b64.txt (data URI for the page)."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
import base64
x = np.linspace(-3, 1.5, 360); y = np.linspace(-1.9, 1.9, 304)
X, Y = np.meshgrid(x, y); Z = X + 1j * Y
L = 1 + Z
h = (np.angle(L) / (2 * np.pi)) % 1.0
m = np.log(np.abs(L) + 1e-12) / np.log(1.5)
band = 0.72 + 0.28 * (m - np.floor(m))
ang = (np.angle(L) / (2 * np.pi) * 12) % 1.0
ray = 0.9 + 0.1 * np.minimum(ang, 1 - ang) * 2
rgb = hsv_to_rgb(np.dstack([h, 0.62 * np.ones_like(h), band * ray]))
fig = plt.figure(figsize=(3.6, 3.04), dpi=100)
ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(rgb, origin="lower", extent=[x[0], x[-1], y[0], y[-1]])
ax.axis("off")
fig.savefig("figures/phase.png", dpi=100)
from PIL import Image
Image.open("figures/phase.png").convert("RGB").save("figures/phase.jpg", quality=78)
b = base64.b64encode(open("figures/phase.jpg", "rb").read()).decode()
open("figures/phase_b64.txt", "w").write("data:image/jpeg;base64," + b)
print(len(b))
