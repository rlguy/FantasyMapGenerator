#include "config.h"

namespace gen {
namespace config {

unsigned int seed = 0;
double resolution = 0.08;
std::string outfileExt = ".png";
std::string outfile = "output" + outfileExt;
double erosionAmount = -1.0;
int erosionIterations = 3;
int imageWidth = 1920;
int imageHeight = 1080;
double defaultExtentsHeight = 20.0;
double drawScale = 1.0;
bool verbose = false;

void print(std::string msg) {
    if (gen::config::verbose) {
        std::cout << msg << std::endl;
    }
}

bool parseOptions(int argc, char **argv) {
    OptionArgs opts;
    void *argtable[] = {
        opts.help       = arg_litn("h", "help", 0, 1, "display this help and exit"),
        opts.seed       = arg_strn("s", "seed", "<uint>", 0, 1, "set random generator seed"),
        opts.timeseed   = arg_litn(NULL, "timeseed", 0, 1, "set seed from system time"),
        opts.resolution = arg_dbln("r", "resolution", "<float>", 0, 1, "level of map detail"),
        opts.outfile    = arg_filen("o", "output", "filename", 0, 1, "output file"),
        opts.output     = arg_filen(NULL, NULL, "<file>", 0, 1, "output file"),
        opts.eroamount  = arg_dbln("e", "erosion-amount", "<float>", 0, 1, "erosion amount"),
        opts.erosteps   = arg_intn(NULL, "erosion-steps", "<int>", 0, 1, "number of erosion iterations"),
        opts.size       = arg_strn(NULL, "size", "<widthpx:heightpx>", 0, 1, "set output image size"),
        opts.drawscale  = arg_dbln(NULL, "draw-scale", "<float>", 0, 1, "set scale of drawn lines/points"),
        opts.drawinfo   = arg_litn(NULL, "drawing-supported", 0, 1, "display whether drawing is supported and exit"),
        opts.verbose    = arg_litn("v", "verbose", 0, 1, "output additional information to stdout"),
        opts.end        = arg_end(20)
    };

    if (arg_nullcheck(argtable) != 0) {
        std::cout << "error: insufficient memory." << std::endl;
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return false;
    }

    int nerrors = arg_parse(argc,argv,argtable);

    std::string progname("map_generation");
    if (opts.help->count > 0) {
        std::cout << "Usage: " << progname;
        arg_print_syntax(stdout, argtable, "\n");
        std::cout << "\nOptions:\n" << std::endl;
        arg_print_glossary(stdout, argtable, "  %-35s %s\n");
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return false;
    }

    if (nerrors > 0) {
        arg_print_errors(stdout, opts.end, progname.c_str());
        std::cout << "Try '" << progname << " --help' for more information." << std::endl;
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return false;
    }

    if (!_displayInfo(opts)) {
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return false;
    }

    if (!_setOptions(opts)) {
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return false;
    }

    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return true;
}

bool _displayInfo(OptionArgs opts) {
    bool isInfoDisplayed = _displayDrawSupportInfo(opts.drawinfo);

    return !isInfoDisplayed;
}

bool _displayDrawSupportInfo(arg_lit *drawinfo) {
    if (drawinfo->count == 0) {
        return false;
    }

    #ifdef PYTHON_RENDERING_SUPPORTED
        std::string result = "True";
    #else
        std::string result = "False";
    #endif

    std::cout << "--drawing-supported=" << result << std::endl;

    return true;
}

bool _setOptions(OptionArgs opts) {
    if (!_setSeed(opts.timeseed, opts.seed)) { return false; }
    if (!_setResolution(opts.resolution)) { return false; }
    if (!_setOutputFile(opts.outfile, opts.output)) { return false; }
    if (!_setErosionAmount(opts.eroamount)) { return false; }
    if (!_setErosionIterations(opts.erosteps)) { return false; }
    if (!_setImageSize(opts.size)) { return false; }
    if (!_setDrawScale(opts.drawscale)) { return false; }
    if (!_setVerbosity(opts.verbose)) { return false; }

    return true;
}

bool _setSeed(arg_lit *timeseed, arg_str *seed) {
    if (timeseed->count > 0) {
        gen::config::seed = (unsigned int)time(NULL);
    } else if (seed->count > 0) {
        std::istringstream istr(seed->sval[0]);
        istr >> gen::config::seed;
    }

    return true;
}

bool _setResolution(arg_dbl *res) {
    if (res->count == 0) {
        return true;
    }

    double r = res->dval[0];
    if (r <= 0) {
        std::cout << "error: resolution must be greater than zero." << std::endl; 
        std::cout << "resolution: " << r << std::endl;
        return false;
    }

    gen::config::resolution = r;

    return true;
}

bool _setOutputFile(arg_file *outfile1, arg_file *outfile2) {
    if (outfile1->count > 0) {
        gen::config::outfile = outfile1->filename[0];
        gen::config::outfileExt = outfile1->extension[0];
    } else if (outfile2->count > 0) {
        gen::config::outfile = outfile2->filename[0];
        gen::config::outfileExt = outfile2->extension[0];
    }

    return true;
}

bool _setErosionAmount(arg_dbl *amount) {
    if (amount->count == 0) {
        return true;
    }

    double e = amount->dval[0];
    if (e <= 0) {
        std::cout << "error: erosion amount must be greater than zero." << std::endl; 
        std::cout << "erosion amount: " << e << std::endl;
        return false;
    }

    gen::config::erosionAmount = e;

    return true;
}

bool _setErosionIterations(arg_int *iterations) {
    if (iterations->count == 0) {
        return true;
    }

    int n = iterations->ival[0];
    if (n < 0) {
        std::cout << "error: erosion iterations must be greater than or equal to zero." << std::endl; 
        std::cout << "erosion iterations: " << n << std::endl;
        return false;
    }

    gen::config::erosionIterations = n;

    return true;
}

bool _setImageSize(arg_str *size) {
    if (size->count == 0) {
        return true;
    }

    std::string sizestr(size->sval[0]);
    std::string delimiter = ":";
    size_t dpos = sizestr.find(delimiter);
    std::string invalidmsg("error: image size must be in form <widthpx:heighpx>.");
    if (dpos == std::string::npos) {
        std::cout << invalidmsg << std::endl;
        return false;
    }

    std::string widthstr = sizestr.substr(0, dpos);
    std::string heightstr = sizestr.substr(dpos + 1, sizestr.size());
    if (widthstr.size() == 0 || heightstr.size() == 0) {
        std::cout << invalidmsg << std::endl;
        return false;
    }

    int widthpx;
    std::istringstream istreamWidth(widthstr);
    istreamWidth >> widthpx;
    if (istreamWidth.fail() || widthpx <= 0) {
        std::cout << "error: invalid image width value." << std::endl;
        std::cout << "width value: " << widthstr << std::endl;
        return false;
    }

    int heightpx;
    std::istringstream istreamHeight(heightstr);
    istreamHeight >> heightpx;
    if (istreamHeight.fail() || heightpx <= 0) {
        std::cout << "error: invalid image height value." << std::endl;
        std::cout << "height value: " << heightstr << std::endl;
        return false;
    }

    gen::config::imageWidth = widthpx;
    gen::config::imageHeight = heightpx;

    return true;
}

bool _setDrawScale(arg_dbl *drawscale) {
    if (drawscale->count == 0) {
        return true;
    }

    double s = drawscale->dval[0];
    if (s <= 0) {
        std::cout << "error: draw scale must be greater than zero." << std::endl; 
        std::cout << "draw size: " << s << std::endl;
        return false;
    }

    gen::config::drawScale = s;

    return true;
}

bool _setVerbosity(arg_lit *verbose) {
    if (verbose->count > 0) {
        gen::config::verbose = true;
    }

    return true;
}

}
}