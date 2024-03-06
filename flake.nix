{
  description = "Utilities and drivers for manipulators (robotic arms)";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/23.11";
    utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, utils, ... }@inputs: {
    overlays.default = (final: prev: {
      pythonPackagesExtensions = prev.pythonPackagesExtensions ++ [
        (python-final: python-prev: {
          kincpp = python-final.callPackage ./nix/pkgs/kincpp {};
        })
      ];
    });
  } // utils.lib.eachSystem [
    "x86_64-linux" "aarch64-linux"
  ] (system:
    let pkgs = import nixpkgs {
          inherit system;
          overlays = [ self.overlays.default ];
          config.allowUnfree = true;
        };
    in {
      packages.default = pkgs.python3Packages.kincpp;
    });
}
