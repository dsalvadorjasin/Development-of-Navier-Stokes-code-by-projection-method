# Render the EnSight Gold output of NS_lid_driven_cavity to a video.
# Usage (headless):  pvbatch make_video.py [EnsightOutput.case] [frames_dir]
# Produces frames_dir/frame_####.png (one per stored timestep); assemble with
#   ffmpeg -framerate 10 -i frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p velocity.mp4
import os
import sys

from paraview.simple import (
    Calculator, ColorBy, CreateView, EnSightReader, GetAnimationScene,
    GetColorTransferFunction, GetScalarBar, Glyph, Render, SaveScreenshot,
    Show, Text,
)

case_file = sys.argv[1] if len(sys.argv) > 1 else "EnsightOutput.case"
out_dir = sys.argv[2] if len(sys.argv) > 2 else "frames"
os.makedirs(out_dir, exist_ok=True)

reader = EnSightReader(CaseFileName=case_file)
reader.PointArrays = ["velocity", "omega", "press"]
reader.UpdatePipeline()

view = CreateView("RenderView")
view.ViewSize = [1280, 960]
view.UseColorPaletteForBackground = 0
view.Background = [1, 1, 1]
view.OrientationAxesVisibility = 0
view.InteractionMode = "2D"

# Coloured surface by |velocity|
surface = Show(reader, view)
surface.SetRepresentationType("Surface")
ColorBy(surface, ("POINTS", "velocity", "Magnitude"))
lut = GetColorTransferFunction("velocity")
lut.ApplyPreset("Viridis", True)
bar = GetScalarBar(lut, view)
bar.Title = "|velocity|"
bar.ComponentTitle = ""
bar.TitleColor = [0, 0, 0]
bar.LabelColor = [0, 0, 0]
bar.Visibility = 1

# Velocity arrows
glyph = Glyph(Input=reader, GlyphType="2D Glyph")
glyph.OrientationArray = ["POINTS", "velocity"]
glyph.ScaleArray = ["POINTS", "velocity"]
glyph.ScaleFactor = 0.08
glyph.GlyphMode = "All Points"
glyph_disp = Show(glyph, view)
ColorBy(glyph_disp, None)
glyph_disp.DiffuseColor = [0, 0, 0]
glyph_disp.AmbientColor = [0, 0, 0]

# Timestep label
label = Text(Text="")
label_disp = Show(label, view)
label_disp.Color = [0, 0, 0]
label_disp.FontSize = 28
label_disp.WindowLocation = "Upper Left Corner"

view.ResetCamera()
view.CameraParallelProjection = 1

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
timesteps = list(reader.TimestepValues) or [0.0]
print("timesteps:", timesteps)

# Fixed colour range (global max of |velocity|) so frames are comparable
vmax = 0.0
for t in timesteps:
    reader.UpdatePipeline(t)
    vmax = max(vmax, reader.PointData["velocity"].GetRange(-1)[1])
lut.RescaleTransferFunction(0.0, vmax)
glyph.ScaleFactor = 0.0 if vmax == 0.0 else 0.08 / vmax

for i, t in enumerate(timesteps):
    scene.AnimationTime = t
    label.Text = "Lid-driven cavity  -  step %d" % int(t)
    Render(view)
    fname = os.path.join(out_dir, "frame_%04d.png" % i)
    SaveScreenshot(fname, view, ImageResolution=[1280, 960])
    print("wrote", fname)
