from PIL import Image, ImageDraw
import glob

frames = []
imgs = []
for i in range(300,1010,10):
    imgs.append("N_{}.png".format(i))

for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)

frames[0].save('movie.gif', format='GIF', append_images=frames[1:], save_all=True, duration = 280, Loop=0)