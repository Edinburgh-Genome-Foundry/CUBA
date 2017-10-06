Example scenario
----------------


=> The following section in HTML-PUG defines the app's layout.


<template lang="pug">
div
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')

  //- Short description of what you can do in this page
  p.center.
    Hi I am an example scenario code for the frontend !
    Duplicate me to get started !

  //- This widget provides links to tweet the page or give feedback
  web-links(:mailSubject="'[CUBA] Feedback on: ' + infos.title",
            tweetMessage="Pre-written message for your tweets !",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")

  //- This is the form, where all the use input goes ! Notice that the v-model
  //- for each form component is of the form ```form.XXXX``
  .form

    h4.formlabel Provide some files to do whatever
    filesuploader(v-model='form.files', :multiple='true',
                  text="Drop multiple files here (or click to select)")

    h4.formlabel Provide some files to do whatever
    el-input(v-model='inputText')

  //- This widget provides the "Submit" button, and manages the request to the
  //- backend. Super important !
  backend-querier(:form='form',
                  :backendUrl='infos.backendUrl',
                  :validateForm='validateForm',
                  submitButtonText='Predict final construct(s)',
                  v-model='queryStatus')

  //- An error message appears if anything goes wrong.
  el-alert(v-if='queryStatus.requestError && !queryStatus.polling.inProgress',
           :title="queryStatus.requestError",
           type="error",
           :closable="false")

  //- The results section appears when the request returned interesting results.
  .results(v-if='!queryStatus.polling.inProgress')
    //- Download button for the file produced by the backend
    download-button(v-if='queryStatus.result.file',
                    :filedata='queryStatus.result.file')
    .results-summary(v-if='queryStatus.result.preview',
                     v-html="queryStatus.result.preview.html")

  //- This widget advertizes for the backend software used in this particular page.
  powered-by(:softwareNames='infos.poweredby')
</template>


=> The following section in Javascript defines the app's variables and logics.


<script>

//- Infos on this scenario.
var infos = {
  title: 'Example Scenario',
  navbarTitle: 'Example Scenario',
  path: 'example_scenario',
  description: '',
  backendUrl: 'start/example_scenario',
  icon: require('assets/images/assembly.svg'),
  poweredby: ['dnacauldron', 'dnafeaturesviewer']
}

export default {
  data: function () {
    return {
      infos: infos,
      form: {
        files: [],
        textInput: ''
      },
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      }
    }
  },
  infos: infos,
  methods: {
    validateForm: function () {
      var errors = []
      if (this.form.files.length === 0) {
        errors.push('Provide at least 2 files: one or more parts and a receptor.')
      }
      if (this.form.inputText.length === 0) {
        errors.push('Provide some text input.')
      }
      return errors
    }
  }
}
</script>


=> The following section applies custome CSS to this particular app.


<style lang='scss' scoped>

h4.formlabel {
  text-align: center;
  text-transform: uppercase;
  margin-top: 40px
}

.form {
  margin: 50px auto;
  max-width: 500px;
}

.title-img {
  height:80px;
  margin-top: -20px;
  margin-bottom: 20px;
}

</style>
