
def rgb2float(r, g, b): return (round(r/255.0, 4), round(g/255.0, 4), round(b/255.0, 4))

# viridis color palette: https://waldyrious.net/viridis-palette-generator/
viridis_16 = [
	(253, 231, 37),
	(210, 226, 27),
	(165, 219, 54),
	(122, 209, 81),
	(84, 197, 104),
	(53, 183, 121),
	(34, 168, 132),
	(31, 152, 139),
	(35, 136, 142),
	(42, 120, 142),
	(49, 104, 142),
	(57, 86, 140) ,
	(65, 68, 135) ,
	(71, 47, 125) ,
	(72, 26, 108) ,
	(68, 1, 84)
]

viridis_32 = [
	(253, 231, 37),
	(234, 229, 26),
	(213, 226, 26),
	(192, 223, 37),
	(168, 219, 52),
	(147, 215, 65),
	(127, 211, 78),
	(108, 205, 90),
	(88, 199, 101),
	(72, 193, 110),
	(58, 186, 118),
	(46, 179, 124),
	(37, 171, 130),
	(32, 163, 134),
	(30, 156, 137),
	(31, 148, 140),
	(34, 140, 141),
	(37, 132, 142),
	(40, 125, 142),
	(43, 117, 142),
	(46, 109, 142),
	(50, 101, 142),
	(54, 93, 141),
	(58, 84, 140),
	(62, 74, 137),
	(66, 65, 134),
	(69, 56, 130),
	(71, 46, 124),
	(72, 35, 116),
	(72, 24, 106),
	(71, 13, 96),
	(68, 1, 84),
]

if __name__ == "__main__":
	# iterate from dark to bright
	for i, c in enumerate(reversed(viridis_32)):
		fc = rgb2float(*c)
		if i == len(viridis_32)-1:
			print(f"vec3({fc[0]}, {fc[1]}, {fc[2]})\t// {c[0]}, {c[1]}, {c[2]}")
		else:
			print(f"vec3({fc[0]}, {fc[1]}, {fc[2]}),\t// {c[0]}, {c[1]}, {c[2]}")