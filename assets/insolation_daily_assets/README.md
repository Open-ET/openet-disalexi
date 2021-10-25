### Deploying the cloud function

Before deploying or calling the cloud functions, the "project" can be set once with the following call, or passed to each gcloud call.
```
gcloud config set project openet-dri
```

The following are the parameters that were set when deploying the function for the first time.  Subsequent deployments only need the project if not set above.
```
gcloud functions deploy cfsr-insolation-daily --project openet-dri --runtime python37 --entry-point cron_scheduler --trigger-http --allow-unauthenticated --memory 512 --timeout 240 --max-instances 1 --service-account="openet-dri@appspot.gserviceaccount.com"
```

### Calling the cloud function

The functions can be called by passing JSON data to the function.
```
gcloud functions call cfsr-insolation-daily --project openet-dri --data '{"date":"2021-10-01"}'
gcloud functions call cfsr-insolation-daily --project openet-dri --data '{"start":"2021-10-01","end":"2021-10-05"}'
```

If no arguments are passed to the scheduler it will check the last 60 days for missing assets.
```
gcloud functions call cfsr-insolation-daily --project openet-dri
```

### Scheduling the job

```
gcloud scheduler jobs update http cfsr-insolation-daily --schedule "0 8 * * *" --uri "https://us-central1-openet.cloudfunctions.net/cfsr-insolation-daily" --description "Update Daily Insolation" --http-method POST --time-zone "UTC" --project openet-dri --max-retry-attempts 5
```
