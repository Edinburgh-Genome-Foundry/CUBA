<template lang="pug">
el-collapse
  el-collapse-item(:title='title', name='1')
    slot Actually there is nothing to learn more
</template>
<script>
export default {
  props: {
    title: {default: ''}
  }
}
</script>
<style scoped>
.el-collapse {
  border: none;
}

.el-collapse-item__header {
  padding-left: 0px;
}
</style>
