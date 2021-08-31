pipeline {
    agent any
    options {
        timeout(time: 1, unit: 'HOURS') 
    }


    stages {

        stage("Creating dist") {
            when { branch 'master' }
            steps {
                sh "chmod 755 ./cicd/script.sh && ./cicd/script.sh"
                // sh "chmod 755 ./cicd/ssm.py && cd cicd && python3 ./ssm.py"
                // office365ConnectorSend 'https://aidashinc.webhook.office.com/webhookb2/818bed7a-c18e-449f-a1e6-4284fdd57ade@13cee7a3-b15d-47d9-8f30-f5bd649e3967/JenkinsCI/ab0b8f048ac24fdd97649d4993b6fd85/0f7eda48-8fed-4b2d-81ef-707d6c02ae75'

            }
        }

    }
    post {
        always {
            
            echo 'One way or another, I have finished'
            deleteDir() /* clean up our workspace */
        }
    }
}


