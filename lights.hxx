#pragma once

#include <vector>
#include <cmath>
#include "math.hxx"
#include "utils.hxx"

class AbstractLight
{
public:

	virtual Vec3f sampleDir(const Vec3f& aSurfPt, float& oLightDist) const
	{
		return Vec3f(0);
	}

	virtual Vec3f getRadiance() const
	{
		return Vec3f(0);
	}

	virtual float getPDF(Vec3f wig, float dist) const
	{
		return 1;
	}

	bool isPoint;
};

//////////////////////////////////////////////////////////////////////////

unsigned int phash = 0x73a08bf0; 

inline float getRand()	// returns random values from 0 to 1
	{
		unsigned int a = phash;
		a = (a ^ 61) ^ (a >> 16);
		a = a + (a << 3);
		a = a ^ (a >> 4);
		a = a * 0x27d4eb2d;
		a = a ^ (a >> 15);
		phash = a;
		return (1.0f*a)/0xffffffff;
	}


class AreaLight : public AbstractLight
{
public:

    AreaLight(
        const Vec3f &aP0,
        const Vec3f &aP1,
        const Vec3f &aP2)
    {
        p0 = aP0;
        e1 = aP1 - aP0;
        e2 = aP2 - aP0;
		
		isPoint = false;

        Vec3f normal = Cross(e1, e2);
        float len    = normal.Length();
        mInvArea     = 2.f / len;
        mFrame.SetFromZ(normal);
		area = len/2.0f;
    }

	virtual Vec3f sampleDir(const Vec3f& aSurfPt, float& oLightDist) const
	{
		Vec3f p = p0;

		Vec2f r(getRand(), getRand());
		Vec2f rt = SampleUniformTriangle(r);

		p = p0 + rt.x*e1 + rt.y*e2;

		Vec3f dir    = p - aSurfPt;
		if(Dot(dir, mFrame.mZ)>=0) return Vec3f(0);
		oLightDist   = dir.Length();
		
		return dir/oLightDist;
	}

	virtual Vec3f getRadiance() const
	{
		return mRadiance;
	}

	virtual float getPDF(Vec3f wig, float dist) const
	{
		return dist*dist / (area * Dot(-wig, mFrame.mZ));
	}

public:
    Vec3f p0, e1, e2;
    Frame mFrame;
    Vec3f mRadiance;
    float mInvArea;
	float iArea;
	float area;
};

//////////////////////////////////////////////////////////////////////////
class PointLight : public AbstractLight
{
public:

    PointLight(const Vec3f& aPosition)
    {
        mPosition = aPosition;
		isPoint = true;
    }

	virtual Vec3f sampleDir(const Vec3f& aSurfPt, float& oLightDist) const
	{
		Vec3f dir	= mPosition - aSurfPt;
		oLightDist	= dir.Length();
		
		return dir/oLightDist;
	}

	virtual Vec3f getRadiance() const
	{
		return mIntensity;
	}

	virtual float getPDF(Vec3f wig, float dist) const
	{
		return dist*dist;
	}

public:

    Vec3f mPosition;
    Vec3f mIntensity;
};


//////////////////////////////////////////////////////////////////////////
class BackgroundLight : public AbstractLight
{
public:
    BackgroundLight()
    {
        mBackgroundColor = Vec3f(135, 206, 250) / Vec3f(255.f);
    }
	
	virtual Vec3f sampleDir(const Vec3f& aSurfPt, float& oLightDist) const
	{
		Vec2f in = Vec2f(getRand(), getRand());
		Vec3f dir = SampleUniformSphereW(in, 0);
		oLightDist   = 1;
		
		return dir/oLightDist;
	}

	virtual Vec3f getRadiance() const
	{
		return Vec3f(0);//mBackgroundColor;
	}

	virtual float getPDF(Vec3f wig, float dist) const
	{
		return UniformSpherePdfW();
	}
	
public:

    Vec3f mBackgroundColor;
};
