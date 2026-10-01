#!/bin/bash
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Building rust lib"
cd ${BASE_DIR}/wpnet_module/fast_ridge
cargo build --release

echo "Copying lib to target dest"
case "$(uname -s)" in
  Darwin*) LIB_FILE="libfast_ridge.dylib" ;;
  *) LIB_FILE="libfast_ridge.so" ;;
esac

mkdir -p ../gpnet/lib
cp target/release/${LIB_FILE} ../gpnet/lib/fast_ridge.so

echo "Finished"
