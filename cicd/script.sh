source ~/venv/bin/activate
export TWINE_USERNAME=aws
export TWINE_PASSWORD=`aws codeartifact get-authorization-token --domain repository-ds-private --domain-owner 300703960986 --query authorizationToken --output text`
export TWINE_REPOSITORY_URL=`aws codeartifact get-repository-endpoint --domain repository-ds-private --domain-owner 300703960986 --repository aidash-ds-pypi --format pypi --query repositoryEndpoint --output text`



python3 setup.py sdist

twine upload --repository codeartifact dist/*.tar.gz
