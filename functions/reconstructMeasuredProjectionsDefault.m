function rec = reconstructMeasuredProjectionsDefault(projBH, r2Vec, degVec, N1, dt1)
  % reconstructMeasuredProjectionsDefault

  rec = iradon(projBH, degVec, 'linear', 'Hann', 1, N1)/dt1;
end
