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
    int width = 500;
    int height = 500;
    std::vector<glm::vec4> image;
    bool freeze = false;
    int current_frame = 0;
    int current_line = 0;
    int max_frames = 1000;
    int max_bounces = 1;
    float epsilon = 2e-4f;
    glm::mat4 view = glm::mat4(1.0f);
    glm::vec3 ground_color = glm::vec3(1.0f, 1.0f, 1.0f);
    glm::vec3 sky_color = glm::vec3(0.5f, 0.7f, 1.0f);
    bool show_normals = false;

    bool antiAliasingOn = true;
    // Indices to determine displayed material of object
    std::vector<int> sphereMaterials;
    std::vector<int> boxMaterials;
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
