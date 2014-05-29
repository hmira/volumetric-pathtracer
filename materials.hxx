#pragma once 

#include "math.hxx"

class Material
{
public:
    Material()
    {
        Reset();
    }

    void Reset()
    {
        mDiffuseReflectance = Vec3f(0);
        mPhongReflectance   = Vec3f(0);
        mPhongExponent      = 1.f;
    }

	void UpdateP()
	{
		float a = mDiffuseReflectance.Max();
		float b = mPhongReflectance.Max();
		pDiff = (a / (a + b));
		pSpec = (b / (a + b));
		ipDiff = 1/pDiff;
		ipSpec = 1/pSpec;
	}

	Vec3f evalBrdf( const Vec3f& wil, const Vec3f& wol ) const
	{
		return evalBrdfDiff(wil, wol) + evalBrdfSpec(wil, wol);
	}

	Vec3f evalBrdfDiff( const Vec3f& wil, const Vec3f& wol ) const
	{
		if( wil.z <= 0 || wol.z <= 0) return Vec3f(0);

		return mDiffuseReflectance / PI_F;
	}

	Vec3f evalBrdfSpec( const Vec3f& wil, const Vec3f& wol ) const
	{
		if( wil.z <= 0 || wol.z <= 0) return Vec3f(0);

		Vec3f R = ReflectLocal(wil);
		return std::pow( Dot(wol, R), mPhongExponent ) * mPhongReflectance *(mPhongExponent+2) / (2*PI_F);
	}

	float getReflectanceDiff(Vec3f wol ) const
	{
		return mDiffuseReflectance.Max();
	}

	float getReflectanceSpec(Vec3f wol ) const
	{
		return mPhongReflectance.Max();
	}

    Vec3f mDiffuseReflectance;
    Vec3f mPhongReflectance;
    float mPhongExponent;

	float pDiff;
	float pSpec;
	float ipDiff;
	float ipSpec;
};
