/////////////////////////////////////////////
//
// Mesh Simplification Tutorial
//
// (C) by Sven Forstmann in 2014
//
// License : MIT
// http://opensource.org/licenses/MIT
//
//https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification
//
// 5/2016: Chris Rorden created minimal version for OSX/Linux/Windows compile

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <string>
#include <math.h>
#define loopi(start_l,end_l) for ( int i=start_l;i<end_l;++i )
#define loopi(start_l,end_l) for ( int i=start_l;i<end_l;++i )
#define loopj(start_l,end_l) for ( int j=start_l;j<end_l;++j )
#define loopk(start_l,end_l) for ( int k=start_l;k<end_l;++k )

float3 barycentric(const float3 &p, const float3 &a, const float3 &b, const float3 &c){
	float3 v0 = b-a;
	float3 v1 = c-a;
	float3 v2 = p-a;
	float d00 = dot(v0, v0);
	float d01 = dot(v0, v1);
	float d11 = dot(v1, v1);
	float d20 = dot(v2, v0);
	float d21 = dot(v2, v1);
	float denom = d00*d11-d01*d01;
	float v = (d11 * d20 - d01 * d21) / denom;
	float w = (d00 * d21 - d01 * d20) / denom;
	float u = 1.0 - v - w;
	return float3(u,v,w);
}

float3 interpolate(const float3 &p, const float3 &a, const float3 &b, const float3 &c, const float3 attrs[3])
{
	float3 bary = barycentric(p,a,b,c);
	float3 out = float3(0,0,0);
	out = out + attrs[0] * bary.x;
	out = out + attrs[1] * bary.y;
	out = out + attrs[2] * bary.z;
	return out;
}

float min(float v1, float v2) {
	return fmin(v1,v2);
}


class SymetricMatrix {

	public:

	// Constructor

	SymetricMatrix(double c=0) { for(int i = 0; i < 10; i++) m[i] = c;  }
	SymetricMatrix(	double m11, double m12, double m13, double m14,
			            double m22, double m23, double m24,
			                        double m33, double m34,
			                                    double m44) {
			 m[0] = m11;  m[1] = m12;  m[2] = m13;  m[3] = m14;
			              m[4] = m22;  m[5] = m23;  m[6] = m24;
			                           m[7] = m33;  m[8] = m34;
			                                        m[9] = m44;
	}

	// Make plane
	SymetricMatrix(double a,double b,double c,double d)
	{
		m[0] = a*a;  m[1] = a*b;  m[2] = a*c;  m[3] = a*d;
		             m[4] = b*b;  m[5] = b*c;  m[6] = b*d;
		                          m[7 ] =c*c; m[8 ] = c*d;
		                                       m[9 ] = d*d;
	}

	double operator[](int c) const { return m[c]; }

	// Determinant

	double det(	int a11, int a12, int a13,
				int a21, int a22, int a23,
				int a31, int a32, int a33)
	{
		double det =  m[a11]*m[a22]*m[a33] + m[a13]*m[a21]*m[a32] + m[a12]*m[a23]*m[a31]
					- m[a13]*m[a22]*m[a31] - m[a11]*m[a23]*m[a32]- m[a12]*m[a21]*m[a33];
		return det;
	}

	const SymetricMatrix operator+(const SymetricMatrix& n) const
	{
		return SymetricMatrix( m[0]+n[0],   m[1]+n[1],   m[2]+n[2],   m[3]+n[3],
						                    m[4]+n[4],   m[5]+n[5],   m[6]+n[6],
						                                 m[ 7]+n[ 7], m[ 8]+n[8 ],
						                                              m[ 9]+n[9 ]);
	}

	SymetricMatrix& operator+=(const SymetricMatrix& n)
	{
		 m[0]+=n[0];   m[1]+=n[1];   m[2]+=n[2];   m[3]+=n[3];
		 m[4]+=n[4];   m[5]+=n[5];   m[6]+=n[6];   m[7]+=n[7];
		 m[8]+=n[8];   m[9]+=n[9];
		return *this;
	}

	double m[10];
};
///////////////////////////////////////////

namespace Simplify
{
	// Global Variables & Strctures
	enum Attributes {
		NONE,
		NORMAL = 2,
		TEXCOORD = 4,
		COLOR = 8
	};
	struct Triangle { 
		int v[4];
		float err[4];
		int deleted,dirty,attr;
		float3 n;
		float3 uvs[3];
	};
	struct Vertex { 
		float3 p;
		int tstart,tcount;
		SymetricMatrix q;
		int border;
	};
	struct Ref { int tid,tvertex; };
	std::vector<Triangle> triangles;
	std::vector<Vertex> vertices;
	std::vector<Ref> refs;
    std::string mtllib;
    std::vector<std::string> materials;

	// Helper functions

	float vertex_error(SymetricMatrix q, float x, float y, float z);
	float calculate_error(int id_v1, int id_v2, float3 &p_result);
	bool flipped(float3 p,int i0,int i1,Vertex &v0,Vertex &v1,std::vector<int> &deleted);
	void update_uvs(int i0,const Vertex &v,const float3 &p,std::vector<int> &deleted);
	void update_triangles(int i0,Vertex &v,std::vector<int> &deleted,int &deleted_triangles);
	void update_mesh(int iteration);
	void compact_mesh();
	void simplify_mesh(int target_count, float agressiveness, bool verbose);
	obj::TriMesh simplify_TriMesh(struct TriMesh *tri_mesh);
	//
	// Main simplification function
	// target_count  : target nr. of triangles
	// agressiveness : sharpness to increase the threshold.
	//                 5..8 are good numbers
	//                 more iterations yield higher quality
	
    void clear() {
	    triangles.clear();
        vertices.clear();
	    refs.clear();
        materials.clear();
    }

	obj::TriMesh simplify(float agressiveness=7, bool verbose=false) {

        int target_count; 
        if(triangles.size() < 10000) {
            target_count = round((float)triangles.size() * 0.1);
        } else {
            target_count = round((float)triangles.size() * 0.01);
        } 

		printf("Input: %zu vertices, %zu triangles (target %d)\n", Simplify::vertices.size(), Simplify::triangles.size(), target_count);
	
		simplify_mesh(target_count, agressiveness, verbose);
		
		obj::TriMesh simple_mesh;
		simple_mesh.indices.resize(4 * triangles.size());
		simple_mesh.texcoords.resize(vertices.size());

		for(int i = 0; i < triangles.size(); i++) {
			Triangle &t = triangles[i];
			for(int k = 0; k < 4; k++)
				simple_mesh.indices[i * 4 + k] = t.v[k];

			for(int k = 0; k < 3; k++) {
				simple_mesh.texcoords[t.v[k]].x = t.uvs[k].x;
				simple_mesh.texcoords[t.v[k]].y = t.uvs[k].y;
			}
		}
        for (auto& v : vertices) {
			simple_mesh.vertices.emplace_back(v.p);	
        }
		
        // Compute the geometric normals for this mesh
        simple_mesh.face_normals.resize(triangles.size());
        obj::compute_face_normals(simple_mesh.indices, simple_mesh.vertices, simple_mesh.face_normals, 0);
		
        simple_mesh.normals.resize(vertices.size());
        std::fill(simple_mesh.normals.begin(), simple_mesh.normals.end(), float3(0.0f));
        obj::compute_vertex_normals(simple_mesh.indices, simple_mesh.face_normals, simple_mesh.normals, 0);
        
        printf("vertices size%ld simple mesh vertices%ld indices size%ld\n", vertices.size(), simple_mesh.vertices.size(), simple_mesh.indices.size() / 4);
		clear();
        return simple_mesh;
	}

	void simplify_mesh(int target_count, float agressiveness=7, bool verbose=false)
	{
		// init
		loopi(0,triangles.size())
            triangles[i].deleted=0;

		// main iteration loop
		int deleted_triangles=0;
		std::vector<int> deleted0,deleted1;
		int triangle_count=triangles.size();
		//int iteration = 0;
		//loop(iteration,0,100)
		for (int iteration = 0; iteration < 100; iteration ++)
		{
			if(triangle_count-deleted_triangles<=target_count)break;

			// update mesh once in a while
			if(iteration%5==0)
			{
				update_mesh(iteration);
			}

			// clear dirty flag
			loopi(0,triangles.size()) triangles[i].dirty=0;

			//
			// All triangles with edges below the threshold will be removed
			//
			// The following numbers works well for most models.
			// If it does not, try to adjust the 3 parameters
			//
			double threshold = 0.000000001*pow(double(iteration+3),agressiveness);

			// target number of triangles reached ? Then break
			if ((verbose) && (iteration%5==0)) {
				printf("iteration %d - triangles %d threshold %g\n",iteration,triangle_count-deleted_triangles, threshold);
			}

			// remove vertices & mark deleted triangles
			loopi(0,triangles.size())
			{
				Triangle &t=triangles[i];
				if(t.err[3]>threshold) continue;
				if(t.deleted) continue;
				if(t.dirty) continue;

				loopj(0,3)if(t.err[j]<threshold)
				{

					int i0=t.v[ j     ]; Vertex &v0 = vertices[i0];
					int i1=t.v[(j+1)%3]; Vertex &v1 = vertices[i1];
					// Border check
					if(v0.border != v1.border)  continue;
					// Compute vertex to collapse to
					float3 p;
					calculate_error(i0,i1,p);
					deleted0.resize(v0.tcount); // normals temporarily
					deleted1.resize(v1.tcount); // normals temporarily
					// don't remove if flipped
					if( flipped(p,i0,i1,v0,v1,deleted0) ) continue;

					if( flipped(p,i1,i0,v1,v0,deleted1) ) continue;

//					if ( (t.v[3] & TEXCOORD) == TEXCOORD  )
//					{
//						update_uvs(i0,v0,p,deleted0);
//						update_uvs(i0,v1,p,deleted1);
//					}

					// not flipped, so remove edge
					v0.p=p;
					v0.q=v1.q+v0.q;
					//printf("%f %f %f |%lf %lf %lf \n", p.x, p.y, p.z, v0.p.x, v0.p.y, v0.p.z);
					int tstart=refs.size();

					update_triangles(i0,v0,deleted0,deleted_triangles);
					update_triangles(i0,v1,deleted1,deleted_triangles);

					int tcount=refs.size()-tstart;

					if(tcount<=v0.tcount)
					{
						// save ram
						if(tcount)memcpy(&refs[v0.tstart],&refs[tstart],tcount*sizeof(Ref));
					}
					else
						// append
						v0.tstart=tstart;

					v0.tcount=tcount;
					break;
				}
				// done?
				if(triangle_count-deleted_triangles<=target_count)break;
			}
		}
		// clean up mesh
		compact_mesh();
	} //simplify_mesh()

	// Check if a triangle flips when this edge is removed

	bool flipped(float3 p,int i0,int i1,Vertex &v0,Vertex &v1,std::vector<int> &deleted)
	{

		loopk(0,v0.tcount)
		{
			Triangle &t=triangles[refs[v0.tstart+k].tid];
			if(t.deleted)continue;

			int s=refs[v0.tstart+k].tvertex;
			int id1=t.v[(s+1)%3];
			int id2=t.v[(s+2)%3];

			if(id1==i1 || id2==i1) // delete ?
			{
				deleted[k]=1;
				continue;
			}
			float3 d1 = vertices[id1].p-p; d1 = normalize(d1);
			float3 d2 = vertices[id2].p-p; d2 = normalize(d2);
			if(fabs(dot(d1, d2))>0.999) return true;
			float3 n;
			n = cross(d1,d2);
			n = normalize(n);
			deleted[k]=0;
			if(dot(n, t.n)<0.2) return true;
		}
		return false;
	}

    // update_uvs

	void update_uvs(int i0,const Vertex &v,const float3 &p,std::vector<int> &deleted)
	{
		loopk(0,v.tcount)
		{
			Ref &r=refs[v.tstart+k];
			Triangle &t=triangles[r.tid];
			if(t.deleted)continue;
			if(deleted[k])continue;
			float3 p1=vertices[t.v[0]].p;
			float3 p2=vertices[t.v[1]].p;
			float3 p3=vertices[t.v[2]].p;
			t.uvs[r.tvertex] = interpolate(p,p1,p2,p3,t.uvs);
		}
	}

	// Update triangle connections and edge error after a edge is collapsed

	void update_triangles(int i0,Vertex &v,std::vector<int> &deleted,int &deleted_triangles)
	{
		float3 p;
		loopk(0,v.tcount)
		{
			Ref &r=refs[v.tstart+k];
			Triangle &t=triangles[r.tid];
			if(t.deleted)continue;
			if(deleted[k])
			{
				t.deleted=1;
				deleted_triangles++;
				continue;
			}
			t.v[r.tvertex]=i0;
			t.dirty=1;
			t.err[0]=calculate_error(t.v[0],t.v[1],p);
			t.err[1]=calculate_error(t.v[1],t.v[2],p);
			t.err[2]=calculate_error(t.v[2],t.v[0],p);
			t.err[3]=min(t.err[0],min(t.err[1],t.err[2]));
			refs.push_back(r);
		}
	}

	// compact triangles, compute edge error and build reference list

	void update_mesh(int iteration)
	{
		if(iteration>0) // compact triangles
		{
			int dst=0;
			loopi(0,triangles.size())
			if(!triangles[i].deleted)
			{
				triangles[dst++]=triangles[i];
			}
			triangles.resize(dst);
		}
		//
		// Init Quadrics by Plane & Edge Errors
		//
		// required at the beginning ( iteration == 0 )
		// recomputing during the simplification is not required,
		// but mostly improves the result for closed meshes
		//
		if( iteration == 0 )
		{
			loopi(0,vertices.size())
			vertices[i].q=SymetricMatrix(0.0);

			loopi(0,triangles.size())
			{
				Triangle &t=triangles[i];
				float3 n,p[3];
				loopj(0,3) p[j]=vertices[t.v[j]].p;
				n = cross(p[1]-p[0],p[2]-p[0]);
				n = normalize(n);
				t.n=n;
				loopj(0,3) vertices[t.v[j]].q =
					vertices[t.v[j]].q+SymetricMatrix(n.x,n.y,n.z,dot(-n, p[0]));
			}
			loopi(0,triangles.size())
			{
				// Calc Edge Error
				Triangle &t=triangles[i];float3 p;
				loopj(0,3) t.err[j]=calculate_error(t.v[j],t.v[(j+1)%3],p);
				t.err[3]=min(t.err[0],min(t.err[1],t.err[2]));
			}
		}

		// Init Reference ID list
		loopi(0,vertices.size())
		{
			vertices[i].tstart=0;
			vertices[i].tcount=0;
		}
		loopi(0,triangles.size())
		{
			Triangle &t=triangles[i];
			loopj(0,3) vertices[t.v[j]].tcount++;
		}
		int tstart=0;
		loopi(0,vertices.size())
		{
			Vertex &v=vertices[i];
			v.tstart=tstart;
			tstart+=v.tcount;
			v.tcount=0;
		}

		// Write References
		refs.resize(triangles.size()*3);
		loopi(0,triangles.size())
		{
			Triangle &t=triangles[i];
			loopj(0,3)
			{
				Vertex &v=vertices[t.v[j]];
				refs[v.tstart+v.tcount].tid=i;
				refs[v.tstart+v.tcount].tvertex=j;
				v.tcount++;
			}
		}

		// Identify boundary : vertices[].border=0,1
		if( iteration == 0 )
		{
			std::vector<int> vcount,vids;

			loopi(0,vertices.size())
				vertices[i].border=0;

			loopi(0,vertices.size())
			{
				Vertex &v=vertices[i];
				vcount.clear();
				vids.clear();
				loopj(0,v.tcount)
				{
					int k=refs[v.tstart+j].tid;
					Triangle &t=triangles[k];
					loopk(0,3)
					{
						int ofs=0,id=t.v[k];
						while(ofs<vcount.size())
						{
							if(vids[ofs]==id)break;
							ofs++;
						}
						if(ofs==vcount.size())
						{
							vcount.push_back(1);
							vids.push_back(id);
						}
						else
							vcount[ofs]++;
					}
				}
				loopj(0,vcount.size()) if(vcount[j]==1)
					vertices[vids[j]].border=1;
			}
		}
	}

	// Finally compact mesh before exiting

	void compact_mesh()
	{
		int dst=0;
		loopi(0,vertices.size())
		{
			vertices[i].tcount=0;
		}
		loopi(0,triangles.size())
		if(!triangles[i].deleted)
		{
			Triangle &t=triangles[i];
			triangles[dst++]=t;
			loopj(0,3)vertices[t.v[j]].tcount=1;
		}
		triangles.resize(dst);
		dst=0;
		loopi(0,vertices.size())
		if(vertices[i].tcount)
		{
			vertices[i].tstart=dst;
			vertices[dst].p=vertices[i].p;
			dst++;
		}
		loopi(0,triangles.size())
		{
			Triangle &t=triangles[i];
			loopj(0,3)t.v[j]=vertices[t.v[j]].tstart;
		}
		vertices.resize(dst);
	}

	// Error between vertex and Quadric

	float vertex_error(SymetricMatrix q, float x, float y, float z)
	{
 		return   q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x + q[4]*y*y
 		     + 2*q[5]*y*z + 2*q[6]*y + q[7]*z*z + 2*q[8]*z + q[9];
	}

	// Error for one edge

	float calculate_error(int id_v1, int id_v2, float3 &p_result)
	{
		// compute interpolated vertex

		SymetricMatrix q = vertices[id_v1].q + vertices[id_v2].q;
		bool   border = vertices[id_v1].border & vertices[id_v2].border;
		float error=0;
		float det = q.det(0, 1, 2, 1, 4, 5, 2, 5, 7);
		if ( det != 0 && !border )
		{

			// q_delta is invertible
			p_result.x = -1/det*(q.det(1, 2, 3, 4, 5, 6, 5, 7 , 8));	// vx = A41/det(q_delta)
			p_result.y =  1/det*(q.det(0, 2, 3, 1, 5, 6, 2, 7 , 8));	// vy = A42/det(q_delta)
			p_result.z = -1/det*(q.det(0, 1, 3, 1, 4, 6, 2, 5,  8));	// vz = A43/det(q_delta)

			error = vertex_error(q, p_result.x, p_result.y, p_result.z);
		}
		else
		{
			// det = 0 -> try to find best result
			float3 p1=vertices[id_v1].p;
			float3 p2=vertices[id_v2].p;
			float3 p3=(p1+p2)/2;
			float error1 = vertex_error(q, p1.x,p1.y,p1.z);
			float error2 = vertex_error(q, p2.x,p2.y,p2.z);
			float error3 = vertex_error(q, p3.x,p3.y,p3.z);
			error = min(error1, min(error2, error3));
			if (error1 == error) p_result=p1;
			if (error2 == error) p_result=p2;
			if (error3 == error) p_result=p3;
		}
		return error;
	}

	char *trimwhitespace(char *str)
	{
		char *end;

		// Trim leading space
		while(isspace((unsigned char)*str)) str++;

		if(*str == 0)  // All spaces?
		return str;

		// Trim trailing space
		end = str + strlen(str) - 1;
		while(end > str && isspace((unsigned char)*end)) end--;

		// Write new null terminator
		*(end+1) = 0;

		return str;
	}

	//Option : Load OBJ
	bool load_obj(const char* filename, obj::MaterialLib &mtl_lib, bool process_uv=false){
		vertices.clear();
		triangles.clear();
		//printf ( "Loading Objects %s ... \n",filename);
		FILE* fn;
		if(filename==NULL)		return -1;
		if((char)filename[0]==0)	return -1;
		if ((fn = fopen(filename, "rb")) == NULL)
		{
			printf ( "File %s not found!\n" ,filename );
			return 0;
		}
		char line[1000];
		memset ( line,0,1000 );
		//int material = -1;
        int cur_mtl = 0;
		std::map<std::string, int> material_map;
		std::vector<float3> uvs;
		std::vector<std::vector<int> > uvMap;

		while(fgets( line, 1000, fn ) != NULL)
		{
			Vertex v;
			float3 uv;

			if (strncmp(line, "mtllib", 6) == 0)
			{
				mtllib = trimwhitespace(&line[7]);
			}
			if (strncmp(line, "usemtl", 6) == 0)
			{
			//	std::string usemtl = trimwhitespace(&line[7]);
			//	if (material_map.find(usemtl) == material_map.end())
			//	{
			//		material_map[usemtl] = materials.size();
			//		materials.push_back(usemtl);
			//	}
		    //	material = material_map[usemtl];
				std::string mtl_name = trimwhitespace(&line[7]);
                cur_mtl = mtl_lib.ids.find(mtl_name)->second;
			}

			if ( line[0] == 'v' && line[1] == 't' )
			{
				if ( line[2] == ' ' )
				if(sscanf(line,"vt %f %f",
					&uv.x,&uv.y)==2)
				{
					uv.z = 0;
					uvs.push_back(uv);
				} else
				if(sscanf(line,"vt %f %f %f",
					&uv.x,&uv.y,&uv.z)==3)
				{
					uvs.push_back(uv);
				}
			}
			else if ( line[0] == 'v' )
			{
				if ( line[1] == ' ' )
				if(sscanf(line,"v %f %f %f",
					&v.p.x,	&v.p.y,	&v.p.z)==3)
				{
					vertices.push_back(v);
				}
			}
			int integers[9];
			if ( line[0] == 'f' )
			{
				Triangle t;
				bool tri_ok = false;
                bool has_uv = false;

				if(sscanf(line,"f %d %d %d",
					&integers[0],&integers[1],&integers[2])==3)
				{
					tri_ok = true;
				}else
				if(sscanf(line,"f %d// %d// %d//",
					&integers[0],&integers[1],&integers[2])==3)
				{
					tri_ok = true;
				}else
				if(sscanf(line,"f %d//%d %d//%d %d//%d",
					&integers[0],&integers[3],
					&integers[1],&integers[4],
					&integers[2],&integers[5])==6)
				{
					tri_ok = true;
				}else
				if(sscanf(line,"f %d/%d/%d %d/%d/%d %d/%d/%d",
					&integers[0],&integers[6],&integers[3],
					&integers[1],&integers[7],&integers[4],
					&integers[2],&integers[8],&integers[5])==9)
				{
					tri_ok = true;
					has_uv = true;
				}
				else
				{
					printf("unrecognized sequence\n");
					printf("%s\n",line);
					while(1);
				}
				if ( tri_ok )
				{
                    for(int i = 0; i < 3; i ++)     
					    t.v[i] = integers[i]-1;
					t.attr = 0;

					if ( process_uv && has_uv )
					{
						std::vector<int> indices;
						indices.push_back(integers[6]-1);
						indices.push_back(integers[7]-1);
						indices.push_back(integers[8]-1);
						uvMap.push_back(indices);
						t.attr |= TEXCOORD;
					}

					t.v[4] = cur_mtl;
					//geo.triangles.push_back ( tri );
					triangles.push_back(t);
					//state_before = state;
					//state ='f';
				}
			}
		}

		if ( process_uv && uvs.size() )
		{
			loopi(0,triangles.size())
			{
				loopj(0,3)
				triangles[i].uvs[j] = uvs[uvMap[i][j]];
			}
		}

		fclose(fn);
        return 1;
		//printf("load_obj: vertices = %lu, triangles = %lu, uvs = %lu\n", vertices.size(), triangles.size(), uvs.size() );
	} // load_obj()

};
