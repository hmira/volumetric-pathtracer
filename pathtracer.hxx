#pragma once

#include <vector>
#include <cmath>
#include <omp.h>
#include <cassert>
#include "renderer.hxx"
#include "rng.hxx"


int sampling_light = 0;

float sigma_a = 0.1f; // ABSORPTION
float sigma_s = 0.2f; // SCATTERING
//float sigma_a = 0.0f; // ABSORPTION
//float sigma_s = 0.0f; // SCATTERING
float sigma_t = sigma_a + sigma_s; // EXTINCTION

class PathTracer : public AbstractRenderer
{
public:

    PathTracer(
        const Scene& aScene,
        int aSeed = 1234
    ) :
        AbstractRenderer(aScene), mRng(aSeed)
    {}

    virtual void RunIteration(int aIteration)
    {
        const int resX = int(mScene.mCamera.mResolution.x);
        const int resY = int(mScene.mCamera.mResolution.y);

        for(int pixID = 0; pixID < resX * resY; pixID++)
        {
            //////////////////////////////////////////////////////////////////////////
            // Generate ray
            const int x = pixID % resX;
            const int y = pixID / resX;

            const Vec2f sample = Vec2f(float(x), float(y)) + mRng.GetVec2f();

            Ray ray = mScene.mCamera.GenerateRay(sample);
            
			// method for path tracing
			Vec3f rad = getLi(ray);
			// method for direct lighting
			//Vec3f rad = directLight(ray);

			mFramebuffer.AddColor(sample, rad);
        }
        mIterations++;
    }

	float selectDistance(const Vec3f& x, const Vec3f& w)
	{
		float xi = getRand();
		return -log(xi) / sigma_t; 
	}

	float e_pow_sigma_t(const float l)
	{
		return exp(-sigma_t * l);
	}

	float e_pow_sigma_t(const Vec3f& x, const Vec3f& y)
	{
		auto l = (x-y).Length();
		return exp(-sigma_t * l);
	}

	/**
	* Counts light contribution to the point in medium. 
	*/
	Vec3f scatteredRadiance(const Vec3f& x)
	{
		// scattering factor, selectDistance pdf, phase function and the light contribution
		return (sigma_s/sigma_t) * ( INV_PI_F * 0.25) * colorSampleLight2(x);
	}

	/**
	* Method for volumetric path tracing.
	*/
	Vec3f getLi(Ray ray)
	{
		Vec3f thrput = Vec3f(1);
		Vec3f accum  = Vec3f(0);
		
		Isect hit;
		bool hits;

		for(int i = 0; i < 1000; i++)
		{
			if(i == 0)
			{
				// inicialization
				hit.dist = 1e36f;
				ray.tmin = EPS_RAY;
				hits = mScene.Intersect(ray, hit);
				if (hit.lightID >= 0) accum += thrput * e_pow_sigma_t(hit.dist) * ((AreaLight*)mScene.GetLightPtr(hit.lightID))->getRadiance();
			}
			
			//if(!hits) break;

			// the surface point hit by the ray
			Vec3f hitPos = ray.org + ray.dir * hit.dist;
			const Material& mat = mScene.GetMaterial( hit.matID );
			Frame frame; frame.SetFromZ(hit.normal);
			const Vec3f wol = frame.ToLocal(-ray.dir);
			
			Vec3f R = ReflectLocal(wol);
			Frame frameR = Frame();
			frameR.SetFromZ(R);

			Vec3f wil, wig;
			float pdf;

			float rand =  getRand();
			// NOTE: there needs to be <=, if there were only <, we get some black dots in the image
			bool sampleDiff = rand <= mat.pDiff;

			/* INCLUDED IN SCATTERING */
			// randomly select a distance on the ray
			auto s = selectDistance(ray.org, ray.dir);
			// if we got shorter distance than ray length, we're in a medium
			if (s < hit.dist)
			{
				// computes the light contribution for the point in medium
				auto xs = ray.org + ray.dir * s;
				accum += thrput * scatteredRadiance(xs);
				wil = sampleMedium(xs, wol, hit, ray, pdf, hits);
				thrput *= (sigma_s / sigma_t);
			}
			// else we got greater distance than ray lenght, we'll take the surface point
			else
			{
				if(!hits) break;
				if(wol.z<=0) break;
				
				// computes the light contribution for the point on surface
				float rho = sampleDiff ? mat.getReflectanceDiff(wol)*mat.ipDiff : mat.getReflectanceSpec(wol)*mat.ipSpec;
				accum += thrput * colorSampleLight(hitPos, wol, frame, mat, sampleDiff, R, frameR);
				wil = sampleBRDF(hitPos, wol, frame, mat, sampleDiff, R, frameR, hit, ray, pdf, hits);
				// russian roulette
				if(rand < rho)
				{
					if(wil.z<0) break;
					thrput *= (sampleDiff ? mat.evalBrdfDiff(wil, wol) : mat.evalBrdfSpec(wil, wol))* wil.z / (rho * pdf);
					//thrput *= mat.evalBrdfDiff(wil, wol) * wil.z / (rho * pdf);
				}
				else 
				{
					// termination of the path
					break;
				}
			}
		}
		return accum;
	}

	/**
	* Method for direct lighting in volumetric renderer.
	*/
	Vec3f directLight(Ray ray)	// poèítá radianci pro smìr, 
	{
		Vec3f accum = Vec3f(0);
		Isect hit;
		hit.dist = 1e36f;
		bool hits = mScene.Intersect(ray, hit);
		if(!hits) return accum;

		Vec3f hitPos = ray.org + ray.dir * hit.dist;
		const Material& mat = mScene.GetMaterial( hit.matID );
		Frame frame; frame.SetFromZ(hit.normal);
		const Vec3f wol = frame.ToLocal(-ray.dir);

		Vec3f R = ReflectLocal(wol);
		Frame frameR = Frame();
		frameR.SetFromZ(R);

		bool sampleDiff = getRand() < mat.pDiff;

		/*INCLUDED IN SCATTERING*/
		// randomly select a distance on the ray
		auto s = selectDistance(ray.org, ray.dir);
		// if we sampled the distance shorter than ray
		if (s < hit.dist)
		{
			// sample the light for the point in medium
			auto xs = ray.org + ray.dir * s;	// point in medium
			accum += scatteredRadiance(xs);
		}
		else
		{
			// we hit light
			if (hit.lightID >= 0) accum += mScene.GetLightPtr(hit.lightID)->getRadiance(); 

			// sample the light for the point on surface
			accum += colorSampleLight(hitPos, wol, frame, mat, sampleDiff, R, frameR);

			Vec3f wil;
			float pdf;
		}

		return accum;
	}

	/**
	* Samples the light for a point in medium.
	*/
	Vec3f colorSampleLight2(Vec3f hitPos)
	{
		Vec3f accum = Vec3f(0);
		const AbstractLight* light =  mScene.GetLightPtr(floor(getRand()*mScene.GetLightCount()));
		float lightDist;
		Vec3f wig = light->sampleDir(hitPos, lightDist);
		if(wig.IsZero()) return accum;

		float p0 = light->getPDF(wig, lightDist) / mScene.GetLightCount();
		if( ! mScene.Occluded(hitPos, wig, lightDist) )
		{
			accum =  light->getRadiance() / p0;
			accum *= e_pow_sigma_t(lightDist);
		}

		return accum;
	}

	/**
	* Samples the light for a surface point.
	*/
	Vec3f colorSampleLight(Vec3f hitPos, Vec3f wol, Frame frame, Material mat, bool sampleDiff, Vec3f R, Frame frameR)
	{
		Vec3f accum = Vec3f(0);

		const AbstractLight* light =  mScene.GetLightPtr(floor(getRand()*mScene.GetLightCount()));
		float lightDist;
		Vec3f wig = light->sampleDir(hitPos, lightDist);
		if(wig.IsZero()) return accum;

		Vec3f wil = frame.ToLocal(wig);
	
		if( ! mScene.Occluded(hitPos, wig, lightDist) )
		{
			float p0 = light->getPDF(wig, lightDist) / mScene.GetLightCount();
			float w = 1.f / p0;
			accum =  w * light->getRadiance() * wil.z * (sampleDiff ? mat.evalBrdfDiff(wil, wol)*mat.ipDiff : mat.evalBrdfSpec(wil, wol)*mat.ipSpec);
			accum *= e_pow_sigma_t(lightDist);
		}

		return accum;
	}

	/**
	* Uniform sampling of the new direction for the point in medium. 
	*/
	Vec3f sampleMedium(
			const Vec3f hitPos, 
			const Vec3f wol,
			Isect & isect,
			Ray & n_ray,
			float & pdf,
			bool & hits
			)
	{
		float p1, w;
		Vec2f in = Vec2f(getRand(), getRand());
		Vec3f wil = SampleUniformSphereW(in, &p1);
		// no need of frame transformation
		Vec3f wig = wil;

		// new ray
		n_ray = Ray(hitPos, wig, 0.00001);
		isect.dist = 1e36f;
		
		// intersecting the scene with the new ray
		hits = mScene.Intersect(n_ray, isect);
		if (hits)
		{
			// if we hit the light
			if (isect.lightID >= 0)
			{	
				const AbstractLight* direct_light = mScene.GetLightPtr(isect.lightID);
			}
		}
		auto wilOut = wil;
		pdf = p1;
		return wilOut;
	}

		
	/**
	* Method for BRDF Importance Sampling. Returns the new direction.
	*/
	Vec3f sampleBRDF(
			const Vec3f hitPos, 
			const Vec3f wol,
			const Frame frame,
			const Material mat,
			const  bool sampleDiff,
			const Vec3f R,
			const Frame frameR,
			Isect & isect,
			Ray & n_ray,
			float & pdf,
			bool & hits
			)
	{
		float p1, w;
		Vec2f in = Vec2f(getRand(), getRand());
		Vec3f wil = sampleDiff ? SampleCosHemisphereW(in, &p1) : SamplePowerCosHemisphereW(in, mat.mPhongExponent, &p1);
		if(!sampleDiff)	wil = frameR.ToWorld(wil);
		Vec3f wig = frame.ToWorld(wil);

		n_ray = Ray(hitPos, wig, 0.00001);
		isect.dist = 1e36f;
		
		hits = mScene.Intersect(n_ray, isect);
		if (hits)
		{
			if (isect.lightID >= 0)
			{	
				const AbstractLight* direct_light = mScene.GetLightPtr(isect.lightID);
			}
		}
		auto wilOut = wil;
		pdf = p1 * (sampleDiff ? mat.pDiff : mat.pSpec);
		return wilOut;
	}





	Vec3f colorSampleBRDF(Vec3f hitPos, Vec3f wol, Frame frame, Material mat, bool sampleDiff, Vec3f R, Frame frameR, Isect & isect, Ray & n_ray, Vec3f & wilOut, float & pdf, bool & hits)
	{
		Vec3f accum = Vec3f(0);

		float p0, p1, w;
		Vec2f in = Vec2f(getRand(), getRand());
		Vec3f wil = sampleDiff ? SampleCosHemisphereW(in, &p1) : SamplePowerCosHemisphereW(in, mat.mPhongExponent, &p1);
		if(!sampleDiff)	wil = frameR.ToWorld(wil);
		Vec3f wig = frame.ToWorld(wil);

		n_ray = Ray(hitPos, wig, 0.00001);
		//Isect isect;
		isect.dist = 1e36f;
		
		hits = mScene.Intersect(n_ray, isect);
		if (hits)
		{
			if (isect.lightID >= 0)
			{	
				const AbstractLight* direct_light = mScene.GetLightPtr(isect.lightID);
				if(!direct_light->isPoint)
				{
					const AreaLight * lght = (AreaLight *) direct_light;
					if(Dot(lght->mFrame.mZ, wig)>0) return accum;
				}
				p0 = direct_light->getPDF(wig, isect.dist) / mScene.GetLightCount();
				//TODO because of eliminating the MIS
				//w = 1 / (p0 + p1);
				w = 1.f / p1;
				accum = w * direct_light->getRadiance() * wil.z * (sampleDiff ? mat.evalBrdfDiff(wil, wol)*mat.ipDiff : mat.evalBrdfSpec(wil, wol)*mat.ipSpec );
			}
		}
		else
		{
			const BackgroundLight* bl = mScene.GetBackground();
			if (bl)
			{
				Vec3f illum =  bl->mBackgroundColor * wil.z;
				Vec3f brdf = (sampleDiff ? mat.evalBrdfDiff(wil, wol) : mat.evalBrdfSpec(wil, wol) );
				accum += ( illum * brdf ) / (p1 * (sampleDiff?mat.pDiff:mat.pSpec));
			}
		}

		wilOut = wil;
		pdf = p1 * (sampleDiff ? mat.pDiff : mat.pSpec);

		return accum;
	}


    Rng mRng;
};
