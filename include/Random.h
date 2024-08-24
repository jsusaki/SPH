#pragma once
#include <random>
#include <chrono>
#include <optional>

#include "Common.h"

class random
{
public:
	random(size_t seed = std::random_device{}()) : mt(seed) {}

public:
	void seed(std::optional<size_t> seed = std::nullopt)
	{
		if (seed.has_value()) 
		{
			mt.seed(seed.value());
		}
		else 
		{
			auto now = std::chrono::high_resolution_clock::now();
			mt.seed(std::random_device{}() ^ now.time_since_epoch().count());
		}
	}

	inline f32 uniform(f32 min, f32 max)
	{
		std::uniform_real_distribution<f32>dist(min, max);
		return dist(mt);
	}

	inline s32 uniformi(s32 min, s32 max)
	{
		std::uniform_int_distribution<s32>dist(min, max);
		return dist(mt);
	}

	inline f32 normal(f32 mu, f32 sigma)
	{
		std::normal_distribution<f32>dist(mu, sigma);
		return dist(mt);
	};

	inline bool bernoulli(f32 p) 
	{
		std::bernoulli_distribution dist(p);
		return dist(mt);
	}

	inline std::vector<f32> uniform(f32 min, f32 max, u32 size)
	{
		std::vector<f32> res(size);
		for (u32 i = 0; i < res.size(); i++)
			res[i] = uniform(min, max);
		return res;
	}

	inline std::vector<s32> uniformi(s32 min, s32 max, u32 size)
	{
		std::vector<s32> res(size);
		for (u32 i = 0; i < res.size(); i++)
			res[i] = uniformi(min, max);
		return res;
	}

	inline std::vector<f32> normal(f32 mu, f32 sigma, u32 size)
	{
		std::vector<f32> res(size);
		for (u32 i = 0; i < res.size(); i++)
			res[i] = normal(mu, sigma);
		return res;
	}

	inline std::vector<f32> bernoulli(f32 p, u32 size)
	{
		std::vector<f32> res(size);
		for (u32 i = 0; i < res.size(); i++)
			res[i] = bernoulli(p);
		return res;
	}

private:
	std::mt19937 mt;
};