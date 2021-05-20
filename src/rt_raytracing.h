#pragma once

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/random.hpp>
#include <cstdlib>
#include <iostream>


#include <vector>
#include <random>
#include <memory>


namespace rt {

struct RTContext {
    int width = 700;
    int height = 700;
    std::vector<glm::vec4> image;
    bool freeze = false;
    int current_frame = 0;
    int current_line = 0;
    int max_frames = 1000;
    int max_bounces = 3;
    float epsilon = 2e-4f;
    glm::mat4 view = glm::mat4(1.0f);
    glm::vec3 ground_color = glm::vec3(1.0f, 0.8f, 1.0f);
    glm::vec3 sky_color = glm::vec3(0.5f, 0.7f, 1.0f);
    bool show_normals = false;
    bool antiAliasingOn = true;

    glm::vec3 groundSphereColor = glm::vec3(0.8, 0.8, 0.0);
    int groundMaterial = 2;
    glm::vec3 lambertianColor = glm::vec3(0.1, 0.2, 0.5);
    float metalFuzz = 0.0;
    glm::vec3 metalColor1 = glm::vec3(0.8, 0.6, 0.2);
    glm::vec3 metalColor2 = glm::vec3(0.77, 0.78, 0.82);
    float dielectricRefraction = 1.5;

    // Indices to determine displayed material of object
    bool showSpheres = true;
    std::vector<int> sphereMaterials;
    bool showBoxes = false;
    std::vector<int> boxMaterials;
    bool showMesh = false;
    int meshMaterial;

    enum materialIndex {LAMBERTIAN, YELLOWMETAL, GREYMETAL, GLASS, HOLLOWGLASS};
};

void setupScene(RTContext &rtx, const char *mesh_filename);
void updateImage(RTContext &rtx, const char* mesh_filename);
void resetImage(RTContext &rtx);
void resetAccumulation(RTContext &rtx);

double random_double();
glm::vec3 random_in_unit_sphere();
glm::vec3 random_in_hemisphere(const glm::vec3& normal);
bool near_zero(glm::vec3 e);
glm::vec3 refract(const glm::vec3& uv, const glm::vec3& n, double etai_over_etat);

}  // namespace rt
