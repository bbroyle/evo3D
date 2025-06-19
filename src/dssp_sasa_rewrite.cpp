/* 
DSSP-derived logic adapted under BSD 2-Clause License
Source: https://github.com/PDB-REDO/dssp

BSD-2-Clause License

Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

// Constants from DSSP
const double PI = 3.141592653589793238462643383279502884;

// Atom radii (in Ångstrom)
const float RADIUS_N = 1.65f;
const float RADIUS_CA = 1.87f;
const float RADIUS_C = 1.76f;
const float RADIUS_O = 1.4f;
const float RADIUS_SIDE_ATOM = 1.8f;
const float RADIUS_WATER = 1.4f;

// Basic 3D point structure
struct Point {
    float x, y, z;
    
    Point() : x(0), y(0), z(0) {}
    Point(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
};

// Vector operations
inline Point operator-(const Point &lhs, const Point &rhs) {
    return Point(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

inline Point operator*(const Point &lhs, float rhs) {
    return Point(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
}

inline float distance_sq(const Point &a, const Point &b) {
    return (a.x - b.x) * (a.x - b.x) +
           (a.y - b.y) * (a.y - b.y) +
           (a.z - b.z) * (a.z - b.z);
}

inline float distance(const Point &a, const Point &b) {
    return std::sqrt(distance_sq(a, b));
}

// Structure to represent a residue
struct Residue {
    int index;
    Point n, ca, c, o;
    std::vector<Point> sideChain;
    float accessibility = 0.0f;
    Point box[2] = {};  // Bounding box: box[0] = min, box[1] = max
    float radius = 0.0f;
    Point center;
    
    Residue(int idx) : index(idx) {
        // Initialize the bounding box
        box[0].x = box[0].y = box[0].z = std::numeric_limits<float>::max();
        box[1].x = box[1].y = box[1].z = -std::numeric_limits<float>::max();
    }
    
    void extendBox(const Point &atom, float radius) {
        if (box[0].x > atom.x - radius) box[0].x = atom.x - radius;
        if (box[0].y > atom.y - radius) box[0].y = atom.y - radius;
        if (box[0].z > atom.z - radius) box[0].z = atom.z - radius;
        if (box[1].x < atom.x + radius) box[1].x = atom.x + radius;
        if (box[1].y < atom.y + radius) box[1].y = atom.y + radius;
        if (box[1].z < atom.z + radius) box[1].z = atom.z + radius;
    }
    
    void finish() {
        // Calculate radius and center of the bounding box
        radius = box[1].x - box[0].x;
        if (radius < box[1].y - box[0].y)
            radius = box[1].y - box[0].y;
        if (radius < box[1].z - box[0].z)
            radius = box[1].z - box[0].z;
            
        center.x = (box[0].x + box[1].x) / 2;
        center.y = (box[0].y + box[1].y) / 2;
        center.z = (box[0].z + box[1].z) / 2;
    }
    
    bool atomIntersectsBox(const Point &atom, float inRadius) const {
        return atom.x + inRadius >= box[0].x &&
               atom.x - inRadius <= box[1].x &&
               atom.y + inRadius >= box[0].y &&
               atom.y - inRadius <= box[1].y &&
               atom.z + inRadius >= box[0].z &&
               atom.z - inRadius <= box[1].z;
    }
};

// Create accumulator for collecting points during surface calculation
class Accumulator {
public:
    struct Candidate {
        Point location;
        double radius;
        double distance;
        
        bool operator<(const Candidate &rhs) const {
            return distance < rhs.distance;
        }
    };
    
    void operator()(const Point &a, const Point &b, double d, double r) {
        double distanceSq = distance_sq(a, b);
        
        d += RADIUS_WATER;
        r += RADIUS_WATER;
        
        double test = d + r;
        test *= test;
        
        if (distanceSq < test && distanceSq > 0.0001) {
            Candidate c = {b - a, r * r, distanceSq};
            m_x.push_back(c);
            std::push_heap(m_x.begin(), m_x.end());
        }
    }
    
    void sort() {
        std::sort_heap(m_x.begin(), m_x.end());
    }
    
    std::vector<Candidate> m_x;
};

// MSurfaceDots class implements a fibonacci sphere for even distribution
class MSurfaceDots {
public:
    static MSurfaceDots &Instance() {
        static MSurfaceDots sInstance;
        return sInstance;
    }
    
    size_t size() const { return points.size(); }
    const Point &operator[](size_t inIx) const { return points[inIx]; }
    double weight() const { return mWeight; }
    
private:
    MSurfaceDots() {
        const int N = 200;  // DSSP (has 200) -- original paper shows 320 (has + or - 1 Ang squared of accuracy)
        auto P = 2 * N + 1;
        const float kGoldenRatio = (1 + std::sqrt(5.0f)) / 2;
        mWeight = (4 * PI) / P;
        
        for (auto i = -N; i <= N; ++i) {
            float lat = std::asin((2.0f * i) / P);
            float lon = static_cast<float>(std::fmod(i, kGoldenRatio) * 2 * PI / kGoldenRatio);
            points.push_back(Point(
                std::sin(lon) * std::cos(lat),
                std::cos(lon) * std::cos(lat),
                std::sin(lat)
            ));
        }
    }
    
    std::vector<Point> points;
    double mWeight;
};

// Function to calculate surface accessibility for a single atom
float calculateAtomSurface(const Point &atom, float atomRadius,
                           const std::vector<Residue*> &neighbors) {
    Accumulator accumulate;
    
    // Find all nearby atoms that might occlude this atom
    for (auto r : neighbors) {
        if (r->atomIntersectsBox(atom, atomRadius)) {
            accumulate(atom, r->n, atomRadius, RADIUS_N);
            accumulate(atom, r->ca, atomRadius, RADIUS_CA);
            accumulate(atom, r->c, atomRadius, RADIUS_C);
            accumulate(atom, r->o, atomRadius, RADIUS_O);
            
            for (const auto &sideAtom : r->sideChain) {
                accumulate(atom, sideAtom, atomRadius, RADIUS_SIDE_ATOM);
            }
        }
    }
    
    accumulate.sort();
    
    float radius = atomRadius + RADIUS_WATER;
    float surface = 0;
    
    MSurfaceDots &surfaceDots = MSurfaceDots::Instance();
    
    // Check each point on the sphere - if not occluded, add to surface area
    for (size_t i = 0; i < surfaceDots.size(); ++i) {
        Point xx = surfaceDots[i] * radius;
        
        bool free = true;
        for (size_t k = 0; free && k < accumulate.m_x.size(); ++k) {
            free = accumulate.m_x[k].radius < distance_sq(xx, accumulate.m_x[k].location);
        }
        
        if (free) {
            surface += static_cast<float>(surfaceDots.weight());
        }
    }
    
    return surface * radius * radius;
}

// Function to calculate surface accessibility for a single residue
float calculateResidueSurface(Residue &res, const std::vector<Residue> &allResidues) {
    std::vector<Residue*> neighbors;
    
    // Find neighboring residues
    for (const auto &r : allResidues) {
        Point center = r.center;
        float radius = r.radius;
        
        if (distance_sq(res.center, center) < (res.radius + radius) * (res.radius + radius)) {
            neighbors.push_back(const_cast<Residue*>(&r));
        }
    }
    
    // Calculate surface for each atom type
    res.accessibility = 
        calculateAtomSurface(res.n, RADIUS_N, neighbors) +
        calculateAtomSurface(res.ca, RADIUS_CA, neighbors) +
        calculateAtomSurface(res.c, RADIUS_C, neighbors) +
        calculateAtomSurface(res.o, RADIUS_O, neighbors);
    
    // Add surface from side chain atoms
    for (const auto &atom : res.sideChain) {
        res.accessibility += calculateAtomSurface(atom, RADIUS_SIDE_ATOM, neighbors);
    }
    
    return res.accessibility;
}

// [[Rcpp::export]]
Rcpp::NumericVector calculateDSSPAccessibility(
    Rcpp::DataFrame atoms,
    Rcpp::DataFrame residues
) {
    /**
     * Calculate solvent accessibility per residue using the DSSP algorithm
     * 
     * @param atoms DataFrame with columns:
     *   - residue_index (int): index matching the residues dataframe
     *   - atom_type (String): N, CA, C, O, or sidechain atom name
     *   - x, y, z (double): Cartesian coordinates
     * 
     * @param residues DataFrame with columns:
     *   - residue_index (int): unique residue identifier
     *   - chain_id (String): chain identifier
     *   - residue_name (String): 3-letter code for amino acid
     *   - residue_number (int): PDB residue number
     * 
     * @return NumericVector with accessibility values for each residue
     */
    
    // Extract data from R
    Rcpp::IntegerVector atom_residue_indices = atoms["residue_index"];
    Rcpp::CharacterVector atom_types = atoms["atom_type"];
    Rcpp::NumericVector atom_x = atoms["x"];
    Rcpp::NumericVector atom_y = atoms["y"];
    Rcpp::NumericVector atom_z = atoms["z"];
    
    Rcpp::IntegerVector residue_indices = residues["residue_index"];
    
    // Create residue objects
    std::vector<Residue> allResidues;
    for (int i = 0; i < residue_indices.size(); i++) {
        allResidues.push_back(Residue(residue_indices[i]));
    }
    
    // Populate residues with atoms
    for (int i = 0; i < atom_residue_indices.size(); i++) {
        int res_idx = -1;
        // Find the corresponding residue by index
        for (size_t j = 0; j < allResidues.size(); j++) {
            if (allResidues[j].index == atom_residue_indices[i]) {
                res_idx = j;
                break;
            }
        }
        
        if (res_idx == -1) {
            Rcpp::stop("Invalid residue index in atoms dataframe");
        }
        
        std::string atom_type = Rcpp::as<std::string>(atom_types[i]);
        Point p(
            (float)atom_x[i],
            (float)atom_y[i],
            (float)atom_z[i]
        );
        
        Residue &res = allResidues[res_idx];
        
        // Add atom to appropriate location in residue
        if (atom_type == "N") {
            res.n = p;
            res.extendBox(p, RADIUS_N + 2 * RADIUS_WATER);
        } 
        else if (atom_type == "CA") {
            res.ca = p;
            res.extendBox(p, RADIUS_CA + 2 * RADIUS_WATER);
        }
        else if (atom_type == "C") {
            res.c = p;
            res.extendBox(p, RADIUS_C + 2 * RADIUS_WATER);
        }
        else if (atom_type == "O") {
            res.o = p;
            res.extendBox(p, RADIUS_O + 2 * RADIUS_WATER);
        }
        else {
            // Assume side chain atom
            res.sideChain.push_back(p);
            res.extendBox(p, RADIUS_SIDE_ATOM + 2 * RADIUS_WATER);
        }
    }
    
    // Finish setting up residues
    for (auto &res : allResidues) {
        res.finish();
    }
    
    // Calculate surface accessibility for each residue
    double totalSurface = 0.0;
    for (auto &res : allResidues) {
        totalSurface += calculateResidueSurface(res, allResidues);
    }
    
    // Return results
    Rcpp::NumericVector result(allResidues.size());
    for (size_t i = 0; i < allResidues.size(); i++) {
        result[i] = allResidues[i].accessibility;
    }
    
    return result;
}