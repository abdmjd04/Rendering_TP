#version 410



#define PI 3.1415926535897932384626433832795


uniform float lightIntensity;
uniform bool blinnPhong;
uniform bool test;
uniform float shininess;
uniform float eta;
uniform sampler2D shadowMap;



const vec4 k_a = vec4(0.2, 0.20, 0.20, 0.0);
const vec4 k_d = vec4(0.4, 0.30, 0.20, 0.0);



in vec4 eyeVector;
in vec4 lightVector;
in vec4 vertColor;
in vec4 vertNormal;
in vec4 lightSpace;

out vec4 fragColor;

//Compute Fresnal

float fresnel(float eta, float cos_theta)
{
    float C_i = sqrt(pow(eta,2.0)-(1-pow(cos_theta,2.0)));
    float F_s = pow(abs((cos_theta-C_i)/(cos_theta+C_i)),2.0);
    float F_p = pow(abs((pow(eta,2.0)*cos_theta-C_i)/(pow(eta,2.0)*cos_theta+C_i)),2.0);
    return 0.5*(F_s+F_p);
}


float X(float cos_theta)
{
    if (cos_theta >= 0 && cos_theta <= 1)
        return 1.0;
    else return 0.0;
}

void main( void )
{


     vec4 ambient_lighting = k_a * vertColor * lightIntensity;

     vec4 diffuse_lighting = k_d * vertColor * max(dot(normalize(vertNormal), normalize(lightVector)),0.0) * lightIntensity;

     vec4 halfVector  = normalize(normalize(eyeVector) + normalize(lightVector));

     float cos_theta = (dot(halfVector,lightVector))/(length(halfVector)*length(lightVector));


     vec4 specular_lighting ;

      // Blinn-Phong Model
     if(blinnPhong){

        specular_lighting =  vertColor * pow(max(dot(normalize(vertNormal), halfVector),0.0),shininess) * lightIntensity;
		
     }


      // Cook-Torrance Model
     else{

		 float alpha        = 1/shininess;
         float alpha_square  = pow(alpha,2.0);
         float cos_theta_i   = dot(normalize(vertNormal), normalize(lightVector));
         float cos_theta_o  = dot(normalize(vertNormal), normalize(eyeVector));
         float cos_theta_h    = dot(halfVector, normalize(vertNormal));

         float tan_square   = (1.0-pow(cos_theta_h,2.0))   / pow(cos_theta_h,2.0);

		 float g_i = 2.0 / (1.0 + sqrt(1.0 + pow(alpha,2.0) * ((1.0-pow(cos_theta_i,2.0))  / pow(cos_theta_i,2.0))));
		 float g_o = 2.0 / (1.0 + sqrt(1.0 + pow(alpha,2.0) * ((1.0-pow(cos_theta_o,2.0))  / pow(cos_theta_o,2.0))));


		 float microfacetNormDist = (X(cos_theta_h)/(PI*pow(cos_theta_h,4.0))) * (alpha_square/(pow(alpha_square + tan_square,2.0)));
         specular_lighting = vertColor * ((microfacetNormDist * g_i * g_o)/ (4.0*cos_theta_i*cos_theta_o))  * lightIntensity;


     //Beckhamann
     //float numerator = exp((pow(cos_theta_h,2.0)-1.0)/(alpha_square * pow(cos_theta_h,2.0)));
     //float d_beck = (numerator)/(PI*alpha_square*pow(cos_theta_h,4.0));
     //specular_lighting = vertColor * ((d_beck * g_i * g_o)/ (4.0*cos_theta_i*cos_theta_o))  * lightIntensity;
     }
    
	 specular_lighting *= fresnel(eta,cos_theta);
        fragColor = (ambient_lighting + diffuse_lighting + specular_lighting);
}
