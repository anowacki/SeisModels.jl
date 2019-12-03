using SeisModels, Test

@testset "IO" begin
    @testset "Mineos" begin
        mktempdir() do dir
            file = joinpath(dir, "model")
            write_mineos(file, AK135)
            ak135 = read_mineos(file)
            @test AK135 ≈ ak135
            @test read_mineos(file) == read_mineos(IOBuffer(String(read(file))))
            let io = IOBuffer()
                write_mineos(io, AK135)
                seekstart(io)
                @test read_mineos(io) ≈ AK135
            end
        end
    end

    @testset "Tvel" begin
        mktempdir() do dir
            file = joinpath(dir, "model.tvel")
            write_tvel(file, AK135)
            @test read_tvel(file) ≈ AK135
            @test read_tvel(file) == read_tvel(IOBuffer(String(read(file))))
            let io = IOBuffer()
                write_tvel(io, AK135)
                seekstart(io)
                @test read_tvel(io) ≈ AK135
            end
            write_tvel(file, AK135, comment1="Line 1", comment2="Line 2")
            lines = readlines(file)
            @test lines[1] == "Line 1"
            @test lines[2] == "Line 2"
        end
    end
end
