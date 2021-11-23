abstract type MaterialStructureType end
abstract type ZincBlende2D <: MaterialStructureType end
abstract type ZincBlende001 <: MaterialStructureType end
abstract type ZincBlende111 <: MaterialStructureType end
abstract type Wurtzite0001 <: MaterialStructureType end
Base.String(::Type{ZincBlende001}) = "ZincBlende001"
Base.String(::Type{ZincBlende111}) = "ZincBlende111"
Base.String(::Type{Wurtzite0001}) = "Wurtzite0001"