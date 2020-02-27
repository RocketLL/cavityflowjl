using FFMPEG, Printf

function generate_gif(src, tgt, fps)
    FFMPEG.exe("-f", "image2", "-framerate", @sprintf("%d", fps), "-i", src, 
    "-filter_complex", "scale=320:-1:flags=lanczos [v], [v] palettegen [p], [v][p] paletteuse", tgt, "-y")
end