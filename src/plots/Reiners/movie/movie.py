from PIL import Image, ImageDraw
import glob

frames = []
imgs = []
for i in range(1,159):
    imgs.append("projected_vel_{}.png".format(i))

for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)

frames[0].save('movie.gif', format='GIF', append_images=frames[1:], save_all=True, duration = 80, Loop=0)