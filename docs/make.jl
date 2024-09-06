using Documenter
using StrainedBandstructures

push!(LOAD_PATH,"../src/")

makedocs(
		modules = [StrainedBandstructures],
		sitename = "StrainedBandStructures.jl",
		authors = "Patricio Farrell, Yiannis Hadjimichael, Christian Merdon",
		format = Documenter.HTML(; repolink = "https://github.com/hadjimy/StrainedBandstructures", mathengine = MathJax3()),
		clean = false,
		source  = "src",
		build   = "build",
		checkdocs = :all,
		doctest = true,
		pages = [
			"Home" => "index.md"
			"Data Types" => [
				"materialdata.md",
				"tensors.md"
			]
		],
	)