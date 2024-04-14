#version 430 core

// Define a uniform struct for lights
struct Light {
    // The matrices are used for shadow mapping. You need to fill it according to how we are filling it when building the normal maps
    // (node_render_shadow_mapping.cpp). 
    // Now, they are filled with identity matrix. You need to modify C++ code in node_render_deferred_lighting.cpp.
    // Position and color are filled.
    mat4 light_projection;
    mat4 light_view;
    vec3 position;
    float radius;
    vec3 color; // Just use the same diffuse and specular color.
    int shadow_map_id;
};

layout(binding = 0) buffer lightsBuffer {
Light lights[4];
};

uniform vec2 iResolution;

uniform sampler2D diffuseColorSampler;
uniform sampler2D normalMapSampler; // You should apply normal mapping in rasterize_impl.fs
uniform sampler2D metallicRoughnessSampler;
uniform sampler2DArray shadow_maps;
uniform sampler2D position;

// uniform float alpha;
uniform vec3 camPos;

uniform int light_count;

layout(location = 0) out vec4 Color;

void main() {
vec2 uv = gl_FragCoord.xy / iResolution;

vec3 pos = texture2D(position,uv).xyz;
vec3 normal = texture2D(normalMapSampler,uv).xyz;

vec4 metalnessRoughness = texture2D(metallicRoughnessSampler,uv);
float metal = metalnessRoughness.x;
float roughness = metalnessRoughness.y;

vec3 result = vec3(0.0);

for(int i = 0; i < light_count; i ++) {

float shadow_map_value = texture(shadow_maps, vec3(uv, lights[i].shadow_map_id)).x;
// HW6_TODO: first comment the line above ("Color +=..."). That's for quick Visualization.
// You should first do the Blinn Phong shading here. You can use roughness to modify alpha. Or you can pass in an alpha value through the uniform above.

Light our_light = lights[i];
// ambient 
float ka = 0.15 / light_count;
vec3 ambient = our_light.color * ka;

//diffusion 
vec3 norm = normalize(normal);
vec3 lightDir = normalize(our_light.position - norm);
float diff = max(dot(norm, lightDir), dot(-norm, lightDir));
float ks = metal * 0.8;
float kd = 1 - ks;
vec3 diffuse = (diff * kd) * our_light.color;

//specular 
vec3 viewDir = normalize(camPos - pos);
vec3 reflecDir = reflect(-lightDir, norm);
float spec = pow(max(dot(viewDir, reflecDir), 0.0), (1 - roughness)*50);
vec3 specular = (spec * our_light.color) * ks;

// shadow 
vec4 posLightSpace = our_light.light_projection * our_light.light_view * vec4(pos, 1.0);
vec3 proj_coord = posLightSpace.xyz / posLightSpace.w;
proj_coord = proj_coord * 0.5 + 0.5;
float closet_depth = texture(shadow_maps, vec3(proj_coord.xy, our_light.shadow_map_id)).r;
float current_depth = posLightSpace.z / posLightSpace.w;
float bias = 0.005;
float shadow = current_depth - bias > closet_depth ? 1.0 : 0.0;

result += ambient +  (1.0 - shadow) * (diffuse + specular);
// After finishing Blinn Phong shading, you can do shadow mapping with the help of the provided shadow_map_value. 
//You will need to refer to the node, node_render_shadow_mapping.cpp, for the light matrices definition. 
//Then you need to fill the mat4 light_projection; mat4 light_view; with similar approach that we fill position and color.
// For shadow mapping, as is discussed in the course, 
//you should compare the value "position depth from the light's view" against the "blocking object's depth.", then you can decide whether it's shadowed.
}
Color = vec4(result, 1.0);
}