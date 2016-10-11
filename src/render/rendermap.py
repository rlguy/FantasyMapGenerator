import cairo
import json
import math

BACKGROUND_RGBA  = (1, 1, 1, 1)
SLOPE_RGBA       = (0, 0, 0, 0.75)
RIVER_RGBA       = (0, 0, 0, 1)
CONTOUR_RGBA     = (0, 0, 0, 1)
BORDER_RGBA      = (0, 0, 0, 1)
CITY_MARKER_RGBA = (0, 0, 0, 1)
TOWN_MARKER_RGBA = (0, 0, 0, 1)
TEXT_RGBA        = (0, 0, 0, 1)

SLOPE_LINE_WIDTH    = 1.0
RIVER_LINE_WIDTH    = 2.5
CONTOUR_LINE_WIDTH  = 1.5
BORDER_LINE_WIDTH   = 6.0
BORDER_DASH_PATTERN = [3, 4]

CITY_MARKER_OUTER_RADIUS = 10
CITY_MARKER_INNER_RADIUS = 5
TOWN_MARKER_RADIUS       = 5

def draw_paths(data, ctx, imgwidth, imgheight):
    for path in data:
        for i in range(0, len(path), 2):
            px = path[i] * imgwidth
            py = imgheight - path[i+1] * imgheight

            if i == 0:
                ctx.move_to(px, py)
            else:
                ctx.line_to(px, py)
    ctx.stroke()

def draw_segments(data, ctx, imgwidth, imgheight):
    for i in range(0, len(data), 4):
        x1 = data[i] * imgwidth
        y1 = imgheight - data[i+1] * imgheight
        x2 = data[i+2] * imgwidth
        y2 = imgheight - data[i+3] * imgheight

        ctx.move_to(x1, y1)
        ctx.line_to(x2, y2)
        ctx.stroke()

def draw_cities(data, ctx, imgwidth, imgheight):
    for i in range(0, len(data), 2):
        px = data[i] * imgwidth
        py = imgheight - data[i+1] * imgheight

        ctx.set_source_rgba(*CITY_MARKER_RGBA)
        ctx.arc(px, py, CITY_MARKER_OUTER_RADIUS, 0, 2*math.pi);
        ctx.fill()

        ctx.set_source_rgba(*BACKGROUND_RGBA)
        ctx.arc(px, py, CITY_MARKER_INNER_RADIUS, 0, 2*math.pi);
        ctx.fill()

def draw_towns(data, ctx, imgwidth, imgheight):
    for i in range(0, len(data), 2):
        px = data[i] * imgwidth
        py = imgheight - data[i+1] * imgheight

        ctx.set_source_rgba(*TOWN_MARKER_RGBA)
        ctx.arc(px, py, TOWN_MARKER_RADIUS, 0, 2*math.pi);
        ctx.fill()

def draw_labels(data, ctx, imgwidth, imgheight):
    for label in data:
        ctx.select_font_face(label["fontface"])
        ctx.set_font_size(int(label["fontsize"]))

        px = label["position"][0] * imgwidth
        py = (1.0 - label["position"][1]) * imgheight

        minx = label["extents"][0] * imgwidth
        miny = (1.0 - label["extents"][1]) * imgheight
        maxx = label["extents"][2] * imgwidth
        maxy = (1.0 - label["extents"][3]) * imgheight
        width = abs(maxx - minx)
        height = abs(maxy - miny)

        ctx.set_source_rgba(*BACKGROUND_RGBA)
        ctx.rectangle(minx, maxy, width, height)
        ctx.fill()

        ctx.set_source_rgba(*TEXT_RGBA)
        ctx.move_to(px, py)
        ctx.show_text(label["text"])

def update_draw_scale(scale):
    global SLOPE_LINE_WIDTH
    global RIVER_LINE_WIDTH
    global CONTOUR_LINE_WIDTH
    global BORDER_LINE_WIDTH
    global BORDER_DASH_PATTERN
    global CITY_MARKER_OUTER_RADIUS
    global CITY_MARKER_INNER_RADIUS
    global TOWN_MARKER_RADIUS

    SLOPE_LINE_WIDTH    *= scale
    RIVER_LINE_WIDTH    *= scale
    CONTOUR_LINE_WIDTH  *= scale
    BORDER_LINE_WIDTH   *= scale
    BORDER_DASH_PATTERN[0] *= scale
    BORDER_DASH_PATTERN[1] *= scale

    CITY_MARKER_OUTER_RADIUS *= scale
    CITY_MARKER_INNER_RADIUS *= scale
    TOWN_MARKER_RADIUS       *= scale

def draw_map(jsonstring, output_filename):
    data = json.loads(jsonstring)

    imgwidth = data["image_width"]
    imgheight = data["image_height"]
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, imgwidth, imgheight)
    ctx = cairo.Context(surface)

    update_draw_scale(data["draw_scale"])
     
    ctx.set_source_rgba(*BACKGROUND_RGBA)
    ctx.rectangle(0, 0, imgwidth, imgheight)
    ctx.fill()

    ctx.set_line_width(SLOPE_LINE_WIDTH)
    ctx.set_source_rgba(*SLOPE_RGBA)
    draw_segments(data["slope"], ctx, imgwidth, imgheight)

    ctx.set_line_width(BORDER_LINE_WIDTH)
    ctx.set_source_rgba(*BACKGROUND_RGBA)
    draw_paths(data["territory"], ctx, imgwidth, imgheight)

    ctx.set_line_width(BORDER_LINE_WIDTH)
    ctx.set_dash(BORDER_DASH_PATTERN)
    ctx.set_line_cap(cairo.LINE_CAP_BUTT)
    ctx.set_line_join(cairo.LINE_JOIN_BEVEL)
    ctx.set_source_rgba(*BORDER_RGBA)
    draw_paths(data["territory"], ctx, imgwidth, imgheight)
    ctx.set_dash([])

    ctx.set_line_width(RIVER_LINE_WIDTH)
    ctx.set_source_rgba(*RIVER_RGBA)
    draw_paths(data["river"], ctx, imgwidth, imgheight)

    ctx.set_line_width(CONTOUR_LINE_WIDTH)
    ctx.set_source_rgba(*CONTOUR_RGBA)
    ctx.set_line_cap(cairo.LINE_CAP_ROUND)
    ctx.set_line_join(cairo.LINE_JOIN_ROUND)
    draw_paths(data["contour"], ctx, imgwidth, imgheight)

    draw_cities(data["city"], ctx, imgwidth, imgheight)
    draw_towns(data["town"], ctx, imgwidth, imgheight)
    draw_labels(data["label"], ctx, imgwidth, imgheight);

    surface.write_to_png(output_filename)
