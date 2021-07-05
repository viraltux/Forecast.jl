using Test
import Forecast: Δt, nΔt

@testset "utils_datetime" begin

    #Date
    ts = [Date("2000-01-01"), Date("2000-01-2")]
    @test_throws AssertionError Δt(ts)

    ts = [Date("2000-01-01"), Date("2000-01-02"),Date("2000-01-03")]
    @test Δt(ts) == Day(1)

    ts = [Date("2000-01-01"), Date("2000-01-11"),Date("2000-01-21")]
    @test Δt(ts) == Day(10)

    ts = [Date("2000-01-04"), Date("2000-01-11"),Date("2000-01-18")]
    @test Δt(ts) == Week(1)    

    ts = [Date("2000-01-01"), Date("2000-02-01"),Date("2000-03-01")]
    @test Δt(ts) == Month(1) 

    ts = [Date("2000-01-01"), Date("2000-03-01"),Date("2000-05-01")]
    @test Δt(ts) == Month(2) 

    ts = [Date("2000-01-01"), Date("2001-01-01"),Date("2002-01-01")]
    @test Δt(ts) == Year(1)
    
    #DateTime
    ts = [DateTime("2000-01-01"), DateTime("2000-01-2")]
    @test_throws AssertionError Δt(ts)


    ts = [DateTime("2000-01-01T00:00:00"),
          DateTime("2000-01-02T00:00:00"),
          DateTime("2000-01-03T00:00:00")]
    @test Δt(ts) == Day(1)

    ts = [DateTime("2000-01-04T00:00:00"),
          DateTime("2000-01-11T00:00:00"),
          DateTime("2000-01-18T00:00:00")]
    @test Δt(ts) == Week(1)
    
    ts = [DateTime("2000-01-01T00:00:00"),
          DateTime("2000-02-01T00:00:00"),
          DateTime("2000-03-01T00:00:00")]
    @test Δt(ts) == Month(1)    

    ts = [DateTime("2000-01-01T00:00:00"),
          DateTime("2001-01-01T00:00:00"),
          DateTime("2002-01-01T00:00:00")]
    @test Δt(ts) == Year(1)    

    ts = [DateTime("2000-01-01T01:00:00"),
          DateTime("2000-01-01T02:00:00"),
          DateTime("2000-01-01T03:00:00")]
    @test Δt(ts) == Hour(1)

    ts = [DateTime("2000-01-01T01:00:00"),
          DateTime("2000-01-01T03:00:00"),
          DateTime("2000-01-01T05:00:00")]
    @test Δt(ts) == Hour(2)

    ts = [DateTime("2000-01-01T00:01:00"),
          DateTime("2000-01-01T00:02:00"),
          DateTime("2000-01-01T00:03:00")]
    @test Δt(ts) == Minute(1)

    ts = [DateTime("2000-01-01T00:00:01"),
          DateTime("2000-01-01T00:00:02"),
          DateTime("2000-01-01T00:00:03")]
    @test Δt(ts) == Second(1)

    ts = [DateTime("2000-01-01T00:00:00.001"),
          DateTime("2000-01-01T00:00:00.002"),
          DateTime("2000-01-01T00:00:00.003")]
    @test Δt(ts) == Millisecond(1)
    
    # Time
    # Time(h, [mi, s, ms, us, ns]) -> Time
    ts = [Time("00:00:00"), Time("01:00:00")]
    @test_throws AssertionError Δt(ts)

    ts = [Time(1,0,0,0,0,0),
          Time(3,0,0,0,0,0),
          Time(5,0,0,0,0,0)]
    @test Δt(ts) == Hour(2)

    ts = [Time(1,0,0,0,0,0),
          Time(2,0,0,0,0,0),
          Time(3,0,0,0,0,0)]
    @test Δt(ts) == Hour(1)

    ts = [Time(0,1,0,0,0,0),
          Time(0,2,0,0,0,0),
          Time(0,3,0,0,0,0)]
    @test Δt(ts) == Minute(1)

    ts = [Time(0,0,1,0,0,0),
          Time(0,0,2,0,0,0),
          Time(0,0,3,0,0,0)]
    @test Δt(ts) == Second(1)

    ts = [Time(0,0,0,1,0,0),
          Time(0,0,0,2,0,0),
          Time(0,0,0,3,0,0)]
    @test Δt(ts) == Millisecond(1)

    ts = [Time(0,0,0,0,1,0),
          Time(0,0,0,0,2,0),
          Time(0,0,0,0,3,0)]
    @test Δt(ts) == Microsecond(1)

    ts = [Time(0,0,0,0,0,0),
          Time(0,0,0,0,0,1),
          Time(0,0,0,0,0,2)]
    @test Δt(ts) == Nanosecond(1)

    # Generic Δt
    ts = [1,2,3]
    @test Δt(ts) == 1

    ts = [1.,3,5]
    @test Δt(ts) == 2.0

    # nΔt
    ts = [Date("2000-01-01"), Date("2000-01-02"), Date("2000-01-03")]
    @test nΔt(ts,  2) == reshape([Date("2000-01-04"),  Date("2000-01-05")],:,1)
    @test nΔt(ts, -2) == reshape([Date("1999-12-30")  Date("1999-12-31")],:,1)
    @test nΔt(ts,  0) == reshape(ts,:,1)

    ts = [DateTime("2000-01-01"), DateTime("2000-01-02"), DateTime("2000-01-03")]
    @test nΔt(ts,  2) == reshape([DateTime("2000-01-04"); DateTime("2000-01-05")],:,1)
    @test nΔt(ts, -2) == reshape([DateTime("1999-12-30"); DateTime("1999-12-31")],:,1)
    @test nΔt(ts,  0) == reshape(ts,:,1)

    ts = [Time("00:00:03"), Time("00:00:04"), Time("00:00:05")]
    @test nΔt(ts,  2) == reshape([Time(0, 0, 6), Time(0, 0, 7)],:,1)
    @test nΔt(ts, -2) == reshape([Time(0, 0, 1), Time(0, 0, 2)],:,1)
    @test nΔt(ts,  0) == reshape(ts,:,1)

    ts = [0,1,2]
    @test nΔt(ts,  2) == reshape([3,4],:,1)
    @test nΔt(ts, -2) == reshape([-2,-1],:,1)
    @test nΔt(ts,  0) == reshape(ts,:,1)

end
