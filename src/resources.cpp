#include "resources.h"

namespace gen{
namespace resources {
    std::string getExecutableDirectory() {
        return std::string(RESOURCES_EXECUTABLE_DIRECTORY);
    }

    std::string getFontDataDirectory() {
        return std::string(RESOURCES_FONT_DATA_DIRECTORY);
    }

    std::string getCityDataDirectory() {
        return std::string(RESOURCES_CITY_DATA_DIRECTORY);
    }

    std::string getFontDataResource() {
        return std::string(RESOURCES_FONT_DATA_RESOURCE);
    }

    std::string getCityDataResource() {
        return std::string(RESOURCES_CITY_DATA_RESOURCE);
    }
}
}