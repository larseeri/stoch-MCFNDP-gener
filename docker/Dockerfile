

## Stage 1: Set up build environment based on a standard image
##############################################################

FROM ubuntu:latest AS builder

## Install build tools and dependencies
RUN apt-get update && apt-get install -y \
        build-essential g++ make

## Set working directory for the build		
##WORKDIR /app

## Copy content from host to image
COPY . .

## Build the application with make tool
RUN make clean HIGHS_SOLVER=1
RUN make HIGHS_SOLVER=1
	
## Stage 2: Set up runtime environment
######################################

## Remark:
## scratch or alpine light weight images cannot support HIGHS
## probably due to missing glibc library and even adding this to alpine:
## RUN apk add --no-cache libc6-compat

FROM ubuntu:latest AS runtime

## Set working directory for the runtime
##WORKDIR /app

## Copy the compiled executable from build stage to image
COPY --from=builder /stoch_gen .

## Copy the default parameter file from the build stage
COPY --from=builder /inout/stoch_config.par .

## Copy the libraries that will used dynamically
COPY --from=builder /HIGHS/HiGHS-1.11.0/install/lib /HIGHS/HiGHS-1.11.0/install/lib

## For diagnostic purposes
##RUN pwd
##RUN ls -al /
##RUN ls -al /HIGHS/HiGHS-1.11.0/install/lib/

# Default command to run the binary when the container starts
CMD ["/stoch_gen", "-h"]

#############################################################################
#############################################################################
## Remark:
##
## Build image named 'stoch_gen_image' with either one of these command line
## instructions:
##
## docker build -t stoch_gen_image .
##
## This empties the cache (i.e. starts building completely anew) and
## displays building process (including output of RUN pwd and RUN ls -al):
##
## docker build --no-cache --progress=plain -t stoch_gen_image .
##
#############################################################################
#############################################################################