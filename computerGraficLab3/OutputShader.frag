uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform vec3 u_pos;
//uniform float u_time;

const float MAX_DIST = 99999.0;

//Поворот матрицы
mat2 rot(float a) {
	float s = sin(a);
	float c = cos(a);
	return mat2(c, -s, s, c);
}

/*Camera*/
struct Camera{
	vec3 Position;
	vec3 Direction;
};
void rot(inout Camera camera){
	camera.Direction.zx *= rot(-u_mouse.y);
	camera.Direction.xy *= rot(u_mouse.x);
}
Camera CameraInit(vec3 Pos, vec3 rd){
	Camera camera;
	camera.Position = Pos;
	camera.Direction = rd;
	return camera;
}


/*Sphere*/
struct Sphere{
	vec3 Position;
	float Radius;
	vec3 Color;
};
Sphere SphereInit(vec3 Pos, float r, vec3 color){
	Sphere sphere;
	sphere.Position = Pos;
	sphere.Radius = r;
	sphere.Color = color;
	return sphere;
}
vec2 sphIntersect(Camera camera,Sphere sphere) {
	vec3 ro = camera.Position - sphere.Position;
	vec3 rd = camera.Direction;
	float ra = sphere.Radius;

	float b = dot(ro, rd);
	float c = dot(ro, ro) - ra * ra;
	float h = b * b - c;
	if(h < 0.0) return vec2(-1.0);
	h = sqrt(h);
	return vec2(-b - h, -b + h);
}



/*Box*/
struct Box{
	vec3 Position;
	vec3 Size;
	vec3 Color;
	vec3 Norm;
};
Box BoxInit(vec3 Pos, vec3 size, vec3 color){
	Box box;
	box.Position = Pos;
	box.Size = size;
	box.Color = color;
	return box;
}
vec2 boxIntersection(Camera camera,out Box box)  {
	vec3 ro = camera.Position - box.Position;
	vec3 rd = camera.Direction;

	vec3 m = 1.0 / rd;

	vec3 n = m * ro;
	vec3 k = abs(m) * box.Size;
	vec3 t1 = -n - k;
	vec3 t2 = -n + k;
	float tN = max(max(t1.x, t1.y), t1.z);
	float tF = min(min(t2.x, t2.y), t2.z);
	if(tN > tF || tF < 0.0) return vec2(-1.0);


	box.Norm = -sign(rd) * step(t1.yzx, t1.xyz) * step(t1.zxy, t1.xyz);
	return vec2(tN, tF);
}

/*Plain*/
float plaIntersect(in vec3 ro, in vec3 rd, in vec4 p) {
	return -(dot(ro, p.xyz) + p.w) / dot(rd, p.xyz);
}

/*Triangle*/

float rayTriangleIntersect(vec3 rayOrigin, vec3 rayDirection, vec3 v0, vec3 v1, vec3 v2, out vec3 normal)
{
    const float EPSILON = 0.000001;
	float t;
    
    vec3 edge1 = v1 - v0;
    vec3 edge2 = v2 - v0;
    vec3 h = cross(rayDirection, edge2);
    float a = dot(edge1, h);
    
    if (abs(a) < EPSILON) return -1.0; // Ray is parallel to triangle
    
    float f = 1.0 / a;
    vec3 s = rayOrigin - v0;
    float u = f * dot(s, h);
    
    if (u < 0.0 || u > 1.0) return -1.0;
    
    vec3 q = cross(s, edge1);
    float v = f * dot(rayDirection, q);
    
    if (v < 0.0 || u + v > 1.0) return -1.0;
    
    t = f * dot(edge2, q);
    
    if (t > EPSILON){
		normal = normalize(cross(edge1, edge2));
        
        
		return t;
	}
    
    return -1.0;
}


/*Tetra*/
struct Tetra{
	vec3 Position;
	float Radius;
	vec3 Norm;
	vec3 Color;
};

Tetra tetraInit(vec3 pos, float r, vec3 color){
	Tetra tetra;
	tetra.Position = pos;
	tetra.Radius = r;
	tetra.Color = color;
	return tetra;
}


vec2 tetraIntersection(Camera camera, out Tetra tetra){
	vec3 ro = camera.Position;
	vec3 rd = camera.Direction;
	float r = tetra.Radius;
	vec3 v1 = tetra.Position + vec3(0.0,(r * 0.57), 0.0);
	vec3 v2 = tetra.Position - vec3(r / 2, r / 3.46,0.0);
	vec3 v3 = tetra.Position + vec3(r / 2, -(r / 3.46),0.0);
	vec3 v4 = tetra.Position - vec3(0.0, 0.0,0.8 * r);

	float it, minIt = MAX_DIST;
	vec3 tn;

	it = rayTriangleIntersect(ro,rd,v1, v2, v3, tn);
	if (it > 0.0 && it < minIt){
		minIt = it;
		tetra.Norm = tn;
	}
	it = rayTriangleIntersect(ro,rd,v2, v1, v4,tn);
	if (it > 0.0 && it < minIt){
		minIt = it;
		tetra.Norm = tn;
	}
	it = rayTriangleIntersect(ro,rd,v1, v3, v4,tn);
	if (it > 0.0 && it < minIt){
		minIt = it;
		tetra.Norm = tn;
	}
	it = rayTriangleIntersect(ro,rd,v3, v2, v4,tn);
	if (it > 0.0 && it < minIt)
	{
		minIt = it;
		tetra.Norm = tn;
	}



	if (minIt == MAX_DIST)
		return vec2(-1.0);

	
	return vec2(minIt);
}





/*Sky*/
struct Sky{
	vec3 light;
	vec3 SkyColor;
	vec3 SunColor;
} sky;
void SkyInit(){
	sky.light = normalize(vec3(-0.5, 0.75, -1.0));
	sky.SkyColor = vec3(0.3, 0.6, 1.0);
	sky.SunColor = vec3(0.95, 0.9, 1.0);
}
vec3 getSky(vec3 rd) {
	vec3 sun = sky.SunColor;
	sun *= max(0.0, pow(dot(rd, sky.light), 32.0));
	return clamp(sun + sky.SkyColor, 0.0, 1.0);
}


vec3 castRay(inout Camera camera) {
	vec3 col;
	vec2 minIt = vec2(MAX_DIST);
	vec2 it;
	vec3 n;

	vec3 ro = camera.Position;
	vec3 rd = camera.Direction;

	Tetra tetra = tetraInit(vec3(-3,-1,1),2.0,vec3(1.0, 1.0, 0.0));
	it = tetraIntersection(camera, tetra);
	if(it.x > 0.0 && it.x < minIt.x) {
		minIt = it;
		n = tetra.Norm;
		col = tetra.Color;
	}

	

	

	Sphere sphere = SphereInit(vec3(0.0, -1.0, 0.0),1,vec3(1.0, 0.2, 0.1));
	it = sphIntersect(camera, sphere);
	if(it.x > 0.0 && it.x < minIt.x) {
		minIt = it;
		vec3 itPos = ro + rd * it.x;
		n = itPos - sphere.Position;
		col = sphere.Color;
	}


	Box box = BoxInit(vec3(0.0, 2.0, 0.0),vec3(1.0),vec3(0.4, 0.6, 0.8));
	it = boxIntersection(camera,box);
	if(it.x > 0.0 && it.x < minIt.x) {
		minIt = it;
		n = box.Norm;
		col = box.Color;
	}
	
	

	vec3 planeNormal = vec3(0.0, 0.0, -1.0);
	it = vec2(plaIntersect(ro, rd, vec4(planeNormal, 1.0)));
	if(it.x > 0.0 && it.x < minIt.x) {
		minIt = it;
		n = planeNormal;
		col = vec3(0.5);
	}



	if(minIt.x == MAX_DIST) return vec3(-1.0);

	float diffuse = max(0.0, dot(sky.light, n));
	float specular = max(0.0, pow(dot(reflect(rd, n), sky.light), 32.0));
	col *= mix(diffuse, specular, 0.5);
	ro += rd * (minIt.x - 0.001);
	//rd = reflect(rd, n);


	camera.Position = ro;
	camera.Direction = rd;

	return col;
}

vec3 traceRay(Camera camera) {
	vec3 col = castRay(camera);
	if(col.x == -1.0) return getSky(camera.Direction);
	Camera newCamer = CameraInit(camera.Position, sky.light);
	if(castRay(newCamer).x != -1.0) col *= 0.5;
	return col;
}

void main() {
	SkyInit();
	//Получение координат пиксиля
	vec2 uv = (gl_TexCoord[0].xy - 0.5) * u_resolution / u_resolution.y;

	Camera camera = CameraInit(u_pos,normalize(vec3(1.0, uv)));
	rot(camera);

	vec3 col = traceRay(camera);

	col.r = pow(col.r, 0.45);
	col.g = pow(col.g, 0.45);
	col.b = pow(col.b, 0.45);
	gl_FragColor = vec4(col, 1.0);
}