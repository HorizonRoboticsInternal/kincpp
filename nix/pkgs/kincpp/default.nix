{ lib
, fetchFromGitHub
, buildPythonPackage
, pythonOlder
, cmake
, numpy
, pybind11
, eigen
}:

let setupPy = ./setup.py;

in buildPythonPackage rec {
  pname = "kincpp";
  version = "1.0.0";
  format = "setuptools";  

  disabled = pythonOlder "3.8";
  
  src = ../../..;

  postPatch = ''
    echo "------------------"
    rm kincpp/__init__.py
    rmdir kincpp
    rm ./setup.py
    ln -s ${setupPy} setup.py
  '';

  propagatedBuildInputs = [
    numpy
  ];

  buildInputs = [
    pybind11
    eigen
  ];

  nativeBuildInputs = [
    cmake
  ];

  dontUseCmakeConfigure = true;  

  pythonImportsCheck = [ "kincpp" ];

  meta = with lib; {
    homepage = "https://github.com/HorizonRoboticsInterna";
    description = ''
      Pybind wrapper functions for forward and inverse kinematics written in C++
    '';
    license = licenses.unfree;
    maintainers = with maintainers; [ breakds ];
    platforms = with platforms; linux;
  };
}
