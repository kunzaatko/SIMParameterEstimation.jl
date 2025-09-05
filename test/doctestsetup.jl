using SIMParameterEstimation
using SIMParameterEstimation: SIMPatternEstimation as SIM_PE
using OffsetArrays: OffsetArrays as OAs

using TestImages
filenames = ["moonsurface.tiff"]; # NOTE: This is a fix for failing doctests since on download, there is a print-out <19-12-24> 
testimage.(filenames; download_only=false);

setup_makiemaestro!() = @eval begin
  using MakieMaestro
  using MakieMaestro.Recipes
  image = Recipes.image
  image! = Recipes.image!
  mosaic = Recipes.mosaic
  nothing
end


return nothing
