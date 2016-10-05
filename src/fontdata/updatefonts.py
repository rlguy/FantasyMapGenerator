import os 
import string
import json
import cairo

VALID_CHARS = (" !\"#$%&'()*+,-./'0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ" + 
               "[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
FORMAT_JSON = False

def updatefonts():
	dir_path = os.path.dirname(os.path.realpath(__file__))
	fontlist_path = os.path.join(dir_path, "fontlist.txt")

	f = open(fontlist_path, 'r')
	lines = f.readlines()
	fonts = []
	for line in lines:
		font_size = line.split(",")
		font_size[0] = font_size[0].strip()
		font_size[1] = int(font_size[1].strip())
		fonts.append(tuple(font_size))

	surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, 1, 1)
	ctx = cairo.Context(surface)

	table = {}
	for font in fonts:
		fontname = font[0]
		fontsize = font[1]

		chardata = {}
		ctx.select_font_face(fontname)
		ctx.set_font_size(fontsize)
		for ch in list(VALID_CHARS):
			x, y, width, height, dx, dy = ctx.text_extents(ch)
			chardata[ch] = [x, y, width, height, dx, dy]

			for i, val in enumerate(chardata[ch]):
				if abs(round(val) - val) < 1e-6:
					chardata[ch][i] = int(round(val))

		if not fontname in table:
			table[fontname] = {}
		table[fontname][str(fontsize)] = chardata

	json_data = json.dumps(table)
	write_path = os.path.join(dir_path, "fontdata.json")
	with open(write_path, 'w') as outfile:
		if FORMAT_JSON:
			json.dump(table, outfile, sort_keys=True, indent=4)
		else:
			json.dump(table, outfile, separators=(',', ':'))
