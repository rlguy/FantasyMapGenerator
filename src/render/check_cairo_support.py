import sys

try:
    import cairo
    is_found = True
except ImportError:
    is_found = False

sys.exit(is_found)
