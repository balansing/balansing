
function [] = mexCompile()
HOME_PATH = '.';
FILE_NAME = 'GenerateEdgesBalanSiNG';
cmd = sprintf('mex(''CXXFLAGS=-ansi -std=c++0x -D_SNU_SOURCE -fPIC -pthread -O3 -DNDEBUG -g'', ''%s/%s.cpp'', ''-largeArrayDims'', ''-outdir'', ''%s'')', HOME_PATH, FILE_NAME, HOME_PATH);
eval(cmd);

HOME_PATH = '.';
FILE_NAME = 'SignedDirectedTriangleEnumeration';
cmd = sprintf('mex(''CXXFLAGS=-ansi -std=c++0x -D_SNU_SOURCE -fPIC -pthread -O3 -DNDEBUG -g'', ''%s/%s.cpp'', ''-largeArrayDims'', ''-outdir'', ''%s'')', HOME_PATH, FILE_NAME, HOME_PATH);
eval(cmd);
end