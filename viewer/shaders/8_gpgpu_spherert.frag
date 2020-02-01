#version 410
#define PI 3.14159265358979323846


uniform mat4 mat_inverse;
uniform mat4 persp_inverse;
uniform sampler2D envMap;
uniform vec3 center;
uniform float radius;


uniform bool bubble;
const int no_bounce = 10;
const float DiameterRatio = 0.99;


uniform float lightIntensity;
uniform bool transparent;
uniform float shininess;
uniform float eta;

in vec4 position;
out vec4 fragColor;

vec3  lightDirection[no_bounce];
vec3  lightNormal[no_bounce];
float lightFresnel[no_bounce];
vec4  lightColour[no_bounce];



vec4 getColorFromEnvironment(in vec3 direction);
bool raySphereIntersect(in vec3 start, in vec3 direction, out vec3 newPoint, float radiusSphere);
float fresnel (vec3 normal, vec3 light, float eta2);
void bubbleTracing(in vec3 eye, in vec3 u, out vec4 colour);
void sphereTracing(in vec3 eye, in vec3 u, out vec4 colour);




vec4 getColorFromEnvironment(in vec3 direction)
{
    float sphereRadius = length(direction);
    float latitute = acos(direction.z/ sphereRadius);
    float longtitute = atan(direction.y , direction.x);
    vec2 coord = vec2(longtitute/(2*PI) + 0.5, latitute/PI);
    return vec4(texture(envMap,coord));
}

float fresnel (vec3 normal, vec3 light, float eta2)
{
    float cosine = dot(normal, light);
    if (cosine < 0)
        eta2 = 1/eta2;
    cosine = abs(cosine);
    float F0 = pow((1 - eta2)/(1 + eta2), 2);
    return (F0 + (1 - F0) * (pow(1 - cosine, 5)));
}


bool raySphereIntersect(in vec3 start, in vec3 direction, out vec3 newPoint, float radiusSphere) {

    //radiusSphere
    vec3 CP = start - center;
    float a = dot(direction, direction);
    float b = 2*dot(direction,CP);
    float c = dot(CP,CP) - (radiusSphere * radiusSphere);
    float delta = pow(b,2) - 4*a*c;


    if (delta > 0){
        float q = (b < 0) ?
                    (-b - sqrt(delta))/2.0:
                    (-b + sqrt(delta))/2.0;

        float t0 = q /a;
        float t1 = c /q;

        // swap if smaller
        if (t0 > t1) {
            float t;
            t = t0;
            t0 = t1;
            t1 = t;
        }

        if (t1 < 0.0)
            return false;

        if (t0 < 0.0)
            newPoint = start + t1*direction;
        else
            newPoint = start + t0*direction;
        return true;
    }
    else
     
        return false;
}

void bubbleTracing(in vec3 eye, in vec3 u, out vec4 colour)
{
    int nextLocation = 2; //1 = env, 2 = glass, 3 = inner air;
    vec3 intersect_point;

    if(raySphereIntersect(eye, u, intersect_point,radius)){

        vec3 vertNormal = normalize(intersect_point - center);
        for(int i = 1; i <no_bounce; i++)
        {
            lightColour[i]  = vec4(0,0,0,0);
            lightFresnel[i] = 0;
        }
        vec3 normal;
        vec3 direction;

        vec3 intersectPointInnerBubble;
        vec3 intersectPointOuter;

        lightFresnel[0]     = fresnel(vertNormal, u ,eta);
        lightNormal[0]      =  vertNormal;
        lightDirection[0]   = normalize(reflect(u, vertNormal));
        lightColour[0]      = getColorFromEnvironment(lightDirection[0]);

        direction = normalize(refract(u, vertNormal, 1/eta));

        //diameter
        float radiusInner = DiameterRatio * radius;

        for(int i = 1; i < no_bounce; i++) {
            normal  = normalize(intersect_point - center);
            if ((dot(direction,normal)<0) && (nextLocation == 2))
            {
                raySphereIntersect(intersect_point+ direction*0.001, direction, intersectPointInnerBubble,radiusInner);
                normal  = normalize(intersectPointInnerBubble - center);
                lightDirection[i+1] = normalize(reflect(direction, normal));
                lightFresnel[i+1]   = fresnel(normal, direction ,eta);
                lightColour[i+1]    = getColorFromEnvironment(lightDirection[i+1]);

                //refract
                lightFresnel[i]     = 1- lightFresnel[i+1] ;
                lightDirection[i]   = normalize(refract(direction, normal, eta));
                lightColour[i]      = getColorFromEnvironment(lightDirection[i]);

                //follow refraction
                direction = lightDirection[i];
                intersect_point = intersectPointInnerBubble;
                nextLocation = 3;
            }
            else if (nextLocation == 3){
                //inner air
                raySphereIntersect(intersect_point+direction*0.001, direction, intersectPointInnerBubble,radiusInner);

                normal              = normalize(center - intersectPointInnerBubble);
                lightDirection[i+1] = normalize(reflect(direction, normal));
                lightFresnel[i+1]   = fresnel(normal, direction ,1/eta);
                lightColour[i+1]    = getColorFromEnvironment(lightDirection[i+1]);

                //refract ouside
                lightFresnel[i]     = 1- lightFresnel[i+1] ;
                lightDirection[i]   = normalize(refract(direction, normal, 1/eta));
                lightColour[i]      = getColorFromEnvironment(lightDirection[i]);

                //follow refraction
                direction           = lightDirection[i];
                intersect_point     = intersectPointInnerBubble;
                nextLocation        = 2;
            }
            else if ((dot(direction,normal)>0) && (nextLocation == 2))
            {
                
                raySphereIntersect(intersect_point + direction*0.001, direction, intersectPointOuter,radius);
                normal              = normalize(center-intersectPointOuter);
                lightDirection[i+1] = normalize(reflect(direction, normal));
                lightFresnel[i+1]   = fresnel(normal, direction ,eta);
                lightColour[i+1]    = getColorFromEnvironment(lightDirection[i+1]);

                //refraction out
                lightFresnel[i]   = 1 - lightFresnel[i+1] ;
                lightDirection[i] = normalize(refract(direction, normal, eta));
                lightColour[i]    = getColorFromEnvironment(lightDirection[i]);

                //follow the reflect
                direction         = lightDirection[i+1];
                intersect_point   = intersectPointOuter;
                nextLocation      = 2;
            }
        }
        for (int i= no_bounce -1; i > 0; i--)
            colour =  lightFresnel[i]*getColorFromEnvironment(lightDirection[i]) + (1-lightFresnel[i]) * colour;
    }
    else
        colour =  getColorFromEnvironment(u) ;
}

void sphereTracing(in vec3 eye, in vec3 u, out vec4 colour)
{

    vec3 intersect_point;

    if(!transparent)
    {
        //metalic effect
        if(raySphereIntersect(eye, u, intersect_point,radius))
        {
            vec3 vertNormal = normalize(intersect_point - center);
            colour =  getColorFromEnvironment(normalize(reflect(u, vertNormal)));
        }
        else
            colour =  getColorFromEnvironment(u) ;
    }
    else
    {
        //refract and reflect
        if(raySphereIntersect(eye, u, intersect_point,radius)){

            vec3 direction;
            vec3 intersect_point2;
            vec3 normal = normalize(intersect_point - center);

            for(int i = 1; i <no_bounce; i++)
            {
                lightColour[i] = vec4(0,0,0,0);
                lightFresnel[i] = 0;
            }

            lightFresnel[0]   = fresnel(normal, u , 1/eta);
            lightNormal[0]    = normal;
            lightDirection[0] = normalize(reflect(u, normal));
            lightColour[0]    = getColorFromEnvironment(lightDirection[0]);

            direction = normalize(refract(u, normal, 1/eta));

            for(int i = 1; i < no_bounce; i++ ){
               
                raySphereIntersect(intersect_point + direction*0.01, direction, intersect_point2,radius);
                normal = normalize(center - intersect_point2);

                //reflect inside 
                lightDirection[i+1] = normalize(reflect(direction, normal));
                lightFresnel[i+1]   = fresnel(normal, direction ,eta);
                lightColour[i+1]    = getColorFromEnvironment(lightDirection[i+1]);

                //refract  ouside 
                lightFresnel[i]   = 1- lightFresnel[i+1] ;
                lightDirection[i] = normalize(refract(direction, normal, eta));
                lightColour[i]    = getColorFromEnvironment(lightDirection[i]);


                //follow reflection
                direction         = lightDirection[i+1];
                intersect_point   = intersect_point2;
            }

            //calculate colour 
            for (int i= no_bounce - 1; i > 0; i--)
                colour = (lightFresnel[i]*getColorFromEnvironment(lightDirection[i])) + (1-lightFresnel[i]) * colour;
        }
        else
            colour =  getColorFromEnvironment(u) ;
    }
}

void main(void)
{
    //pixel coordinates.
    vec4 worldPos = position;
    worldPos.z = 1; // near clipping plane
    worldPos = persp_inverse * worldPos;
    worldPos /= worldPos.w;
    worldPos.w = 0;
    worldPos = normalize(worldPos);
    //keep default eta when start

    // ray direction
    vec3 u = normalize((mat_inverse * worldPos).xyz); // directions
    vec3 eye = (mat_inverse * vec4(0, 0, 0, 1)).xyz; // origin
    vec4 colour = vec4(0,0,0,1);

    if (bubble)
        bubbleTracing(eye,u,colour);
    else
        sphereTracing(eye,u,colour);

    fragColor = colour;
}
