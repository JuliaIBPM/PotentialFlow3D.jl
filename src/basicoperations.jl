cross(r1::SVector{3,T},r2::SVector{3,T}) where {T} = SVector{3,T}([r1[2]*r2[3]-r1[3]*r2[2],
                                                                   r1[3]*r2[1]-r1[1]*r2[3],
                                                                   r1[1]*r2[2]-r1[2]*r2[1]])