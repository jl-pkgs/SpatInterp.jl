@testset "angle2quad" begin
  @test angle2quad(45) == 2
  @test angle2quad(135) == 4
  @test angle2quad(135 + 90 * 1) == 3
  @test angle2quad(135 + 90 * 2) == 1
end
