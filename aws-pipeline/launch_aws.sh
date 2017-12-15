#!/bin/bash

set -e

SUPERBATCHPATH=$1
AWS_ACCESS_KEY=$2
AWS_SECRET_KEY=$3
KEYNAME=$4

usage()
{
    BASE=$(basename -- "$0")
    echo "Launch amazon instance to run MUTEA batch
Usage:
    $BASE <superbatchpath> <aws_access_key> <aws_secret_key> <keyname>
"
    exit 1
}
test -z ${SUPERBATCHPATH} && usage
test -z ${AWS_ACCESS_KEY} && usage
test -z ${AWS_SECRET_KEY} && usage
test -z ${KEYNAME} && usage

# Instance details
SPOT_PRICE=0.05
INSTANCE_TYPE=c4.xlarge
IMAGE_ID=ami-80861296 # TODO change

STARTUP_SCRIPT=$(cat run_from_aws.sh | \
    sed "s/\$1/${AWS_ACCESS_KEY}/" | sed "s~\$2~${AWS_SECRET_KEY}~" | \
    sed "s~\$3~${SUPERBATCHPATH}~")
STARTUP_SCRIPT_ENCODE="$(echo "${STARTUP_SCRIPT}" | base64 -w 0)"

#LAUNCH_SPEC="{\"EbsOptimized\":true, \"ImageId\":\"${IMAGE_ID}\",\"Placement\":{\"AvailabilityZone\": \"us-east-1b\"},\"SecurityGroupIds\":[\"sg-5e914222\"], \"KeyName\":\"${KEYNAME}\",\"InstanceType\":\"${INSTANCE_TYPE}\", \"UserData\":\"${STARTUP_SCRIPT_ENCODE}\", \"BlockDeviceMappings\": [ {\"DeviceName\": \"/dev/sdf\",\"Ebs\": {\"VolumeSize\": 150,\"DeleteOnTermination\": true,\"VolumeType\": \"gp2\"}}]}"

LAUNCH_SPEC="{\"EbsOptimized\":true,\"ImageId\":\"${IMAGE_ID}\",\"Placement\":{\"AvailabilityZone\": \"us-east-1b\"},\"SecurityGroupIds\":[\"sg-5e914222\"], \"KeyName\":\"${KEYNAME}\",\"InstanceType\":\"${INSTANCE_TYPE}\", \"UserData\":\"${STARTUP_SCRIPT_ENCODE}\",\"BlockDeviceMappings\": [ {\"DeviceName\": \"/dev/sdf\",\"Ebs\": {\"VolumeSize\": 150,\"DeleteOnTermination\": true,\"VolumeType\": \"gp2\"}}]}"

aws ec2 request-spot-instances \
    --spot-price ${SPOT_PRICE} \
    --instance-count 1 \
    --type one-time \
    --launch-specification "${LAUNCH_SPEC}"
