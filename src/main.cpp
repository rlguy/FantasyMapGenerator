#include <cmath>

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <random>

#include "mapgenerator.h"
#include "render.h"

double randomDouble(double min, double max) {
    return min + (double)rand() / ((double)RAND_MAX / (max - min));
}

dcel::Point randomPoint(Extents2d &extents) {
    double px = randomDouble(extents.minx, extents.maxx);
    px = randomDouble(extents.minx, extents.maxx);
    double py = randomDouble(extents.miny, extents.maxy);

    return dcel::Point(px, py);
}

dcel::Point randomDirection() {
    double angle = randomDouble(0, 2*3.14159);
    return dcel::Point(sin(angle), cos(angle));
}

void outputMap(gen::MapGenerator &map) {
    std::vector<char> drawdata = map.getDrawData();
    #ifdef PYTHON_RENDERING_SUPPORTED
        std::string filename("output.png");
        gen::render::drawMap(drawdata, filename);

        std::cout << "Wrote map to image: " << filename << std::endl;
    #else
        std::string filename("output.json");
        std::ofstream file(filename);
        file << std::string(drawdata.data());
        file.close();

        std::cout << "Project build without drawing support. " << 
                     "Install Python and the Pycairo graphics library " <<
                     "(http://cairographics.org/pycairo/) to enable drawing " <<
                     "support." << std::endl;
        std::cout << "Wrote map draw data to file: " << filename << std::endl;
    #endif
}

std::vector<std::string> getLabelNames(int num) {
    std::ifstream file(gen::resources::getCityDataResource());
    std::string jsonstr((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());
    jsoncons::json json = jsoncons::json::parse(jsonstr);

    std::vector<std::string> countries;
    for (const auto& member : json.members()) {
        countries.push_back(member.name());
    }

    std::vector<std::string> cities;
    while ((int)cities.size() < num) {
        int randidx = rand() % (int)countries.size();
        std::string country = countries[randidx];
        for (const auto& member : json[country].elements()) {
            cities.push_back(member.as<std::string>());
        }
    }

    std::string temp;
    for (int i = cities.size() - 2; i >= 0; i--) {
        int j = (rand() % (int)(i - 0 + 1));
        temp = cities[i];
        cities[i] = cities[j];
        cities[j] = temp;
    }

    std::vector<std::string> labelnames;
    labelnames.insert(labelnames.end(), cities.begin(), cities.begin() + num);

    return labelnames;
}

int main() {
    time_t seed = time(NULL);
    srand(seed);
    for (int i = 0; i < 1000; i++) { rand(); }
    std::cout << "Generating map with seed value: " << (unsigned int)seed << std::endl;

    Extents2d extents(0, 0, 1.7777*20.0, 20.0);
    gen::MapGenerator map(extents, 0.08);
    map.initialize();
   
    double pad = 5.0;
    Extents2d expandedExtents(extents.minx - pad, extents.miny - pad,
                              extents.maxx + pad, extents.maxy + pad);
    
    int n = randomDouble(100, 250);
    double minr = 1.0;
    double maxr = 8.0;
    for (int i = 0; i < n; i++) {
        dcel::Point p = randomPoint(expandedExtents);
        double r = randomDouble(minr, maxr);
        double strength = randomDouble(0.5, 1.5);
        if (randomDouble(0, 1) > 0.5) {
            map.addHill(p.x, p.y, r, strength);
        } else {
            map.addCone(p.x, p.y, r, strength);
        }
    }

    if (randomDouble(0, 1) > 0.5) {
        dcel::Point p = randomPoint(expandedExtents);
        double r = randomDouble(6.0, 12.0);
        double strength = randomDouble(1.0, 3.0);
        map.addCone(p.x, p.y, r, strength);
    }

    if (randomDouble(0, 1) > 0.1) {
        dcel::Point dir = randomDirection();
        dcel::Point lp = randomPoint(extents);
        double slopewidth = randomDouble(0.5, 5.0);
        double strength = randomDouble(2.0, 3.0);
        map.addSlope(lp.x, lp.y, dir.x, dir.y, slopewidth, strength);
    }
    
    if (randomDouble(0, 1) > 0.5) {
        map.normalize();
    } else {
        map.round();
    }

    if (randomDouble(0, 1) > 0.5) {
        map.relax();
    }
  
    int erosionSteps = 15;
    double erosionAmount = randomDouble(0.2, 0.35);
    for (int i = 0; i < erosionSteps; i++) {
        map.erode(erosionAmount / (double)erosionSteps);
    }
    map.setSeaLevelToMedian();

    int numCities = (int)randomDouble(3, 7);
    int numTowns = (int)randomDouble(8, 25);
    int numLabels = 2*numCities + numTowns;
    std::vector<std::string> labelNames = getLabelNames(numLabels);
    for (int i = 0; i < numCities; i++) {
        std::string cityName = labelNames.back();
        labelNames.pop_back();
        std::string territoryName = labelNames.back();
        labelNames.pop_back();
        std::transform(territoryName.begin(), territoryName.end(), 
                       territoryName.begin(), ::toupper);
        map.addCity(cityName, territoryName);
    }

    for (int i = 0; i < numTowns; i++) {
        std::string townName = labelNames.back();
        labelNames.pop_back();
        map.addTown(townName);
    }

    outputMap(map);
    
    return 0;
}
