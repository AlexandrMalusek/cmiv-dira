function rec = reconstructIteratedProjectionsDefault(proj, r2Vec, degVec, N1, dt1)
  % reconstructMeasuredProjectionsDefault

  rec = iradon(proj, degVec, 'linear', 'Hann', 1, N1)/dt1; 
end
