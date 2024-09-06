using Documenter
using StrainedBandstructures

push!(LOAD_PATH,"../src/")

makedocs(
		modules = [StrainedBandstructures],
		sitename = "StrainedBandStructures.jl",
		authors = "Patricio Farrell, Yiannis Hadjimichael, Christian Merdon",
		repo = "github.com/hadjimy/StrainedBandstructures",
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