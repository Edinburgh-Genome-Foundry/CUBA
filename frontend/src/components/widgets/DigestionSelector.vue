<template lang="pug">
el-select(v-model='value', placeholder='Enzymes in this digestion', filterable, :multiple='true')
  el-option(v-for='enzyme in enzymes.list', :label='enzyme', :value='enzyme', :key='enzyme')
</template>

<script>
import enzymes from '../../assets/enzymes'

export default {
  props: {
    originalValue: {
      default: () => ([])
    }
  },
  data: function () {
    return {
      enzymes: enzymes,
      value: this.originalValue
    }
  },
  watch: {
    value: function (val) {
      this.$emit('input', val)
    }
  }
}
</script>

<style scoped>
.el-tag {
  height: 35px !important;
}

.el-select {
  width: 100%;
}
</style>
