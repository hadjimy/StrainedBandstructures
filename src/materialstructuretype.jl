abstract type MaterialStructureType end
abstract type ZincBlende2D <: MaterialStructureType end
abstract type ZincBlende001 <: MaterialStructureType end
abstract type ZincBlende111_C14 <: MaterialStructureType end
abstract type ZincBlende111_C15 <: MaterialStructureType end
abstract type ZincBlende111_C14_C15 <: MaterialStructureType end
abstract type Wurtzite0001 <: MaterialStructureType end
Base.String(::Type{ZincBlende001}) = "ZincBlende001"
Base.String(::Type{ZincBlende111_C14}) = "ZincBlende111_C14"
Base.String(::Type{ZincBlende111_C15}) = "ZincBlende111_C14"
Base.String(::Type{ZincBlende111_C14_C15}) = "ZincBlende111_C14_C15"
Base.String(::Type{Wurtzite0001}) = "Wurtzite0001"