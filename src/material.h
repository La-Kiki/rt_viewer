#pragma once

#include <glm/glm.hpp>

#include <cstdlib>
#include <iostream>

struct HitRecord;

namespace rt {

class Material {
  public:
     virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const = 0;
};

class Lambertian : public Material {
public:
    Lambertian(const glm::vec3 &a) : albedo(a) {}

    virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const override {
        auto scatter_direction = rec.normal + glm::normalize(rt::random_in_unit_sphere());
        if (rt::near_zero(scatter_direction)) {
            scatter_direction = rec.normal;
        }
            
        scattered = Ray(rec.p, scatter_direction);
        attenuation = albedo;
        return true;
    }

public:
    glm::vec3 albedo;
};

class Metal : public Material {
public:
    Metal(const glm::vec3 &a) : albedo(a) {}

    virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const override {
        glm::vec3 reflected = glm::reflect(glm::normalize(r_in.direction()), rec.normal);
        scattered = Ray(rec.p, reflected);
        attenuation = albedo;
        return (glm::dot(scattered.direction(), rec.normal) > 0);
    }

public:
    glm::vec3 albedo;
};

}  // namespace rt